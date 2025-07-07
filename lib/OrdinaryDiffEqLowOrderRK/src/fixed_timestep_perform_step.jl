function initialize!(integrator, cache::RK4ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RK4ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    halfdt = dt / 2
    k₁ = integrator.fsalfirst
    ttmp = t + halfdt
    k₂ = f(uprev + halfdt * k₁, p, ttmp)
    k₃ = f(uprev + halfdt * k₂, p, ttmp)
    k₄ = f(uprev + dt * k₃, p, t + dt)
    u = uprev + (dt / 6) * (2 * (k₂ + k₃) + (k₁ + k₄))
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    if integrator.opts.adaptive
        # Shampine Solving ODEs and DDEs with Residual Control Estimate
        k₅ = integrator.fsallast

        # one(t) so that types are correct but unitless
        σ₁ = one(t) * (1 // 2) - sqrt(one(t) * 3) / 6
        σ₂ = one(t) * (1 // 2) + sqrt(one(t) * 3) / 6
        p1 = (1 - σ₁) * uprev + σ₁ * u +
             σ₁ * (σ₁ - 1) * ((1 - 2σ₁) * (u - uprev) + (σ₁ - 1) * dt * k₁ + σ₁ * dt * k₅)
        p2 = (1 - σ₂) * uprev + σ₂ * u +
             σ₂ * (σ₂ - 1) * ((1 - 2σ₂) * (u - uprev) + (σ₂ - 1) * dt * k₁ + σ₂ * dt * k₅)
        pprime1 = k₁ +
                  σ₁ * (-4 * dt * k₁ - 2 * dt * k₅ - 6 * uprev +
                   σ₁ * (3 * dt * k₁ + 3 * dt * k₅ + 6 * uprev - 6 * u) + 6 * u) / dt
        pprime2 = k₁ +
                  σ₂ * (-4 * dt * k₁ - 2 * dt * k₅ - 6 * uprev +
                   σ₂ * (3 * dt * k₁ + 3 * dt * k₅ + 6 * uprev - 6 * u) + 6 * u) / dt
        e1 = integrator.opts.internalnorm(
            calculate_residuals(dt * (f(p1, p, t + σ₁ * dt) -
                                      pprime1), uprev, u,
                integrator.opts.abstol,
                integrator.opts.reltol,
                integrator.opts.internalnorm,
                t),
            t)
        e2 = integrator.opts.internalnorm(
            calculate_residuals(dt * (f(p2, p, t + σ₂ * dt) -
                                      pprime2), uprev, u,
                integrator.opts.abstol,
                integrator.opts.reltol,
                integrator.opts.internalnorm,
                t),
            t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
        integrator.EEst = convert(typeof(one(t)), 2.1342) * max(e1, e2)
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

get_fsalfirstlast(cache::RK4Cache, u) = (cache.fsalfirst, cache.k)
function initialize!(integrator, cache::RK4Cache)
    @unpack tmp, fsalfirst, k₂, k₃, k₄, k = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # pre-start FSAL
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RK4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tmp, fsalfirst, k₂, k₃, k₄, k, atmp, stage_limiter!, step_limiter!, thread = cache
    k₁ = fsalfirst
    halfdt = dt / 2
    ttmp = t + halfdt
    @.. broadcast=false thread=thread tmp=uprev + halfdt * k₁
    stage_limiter!(tmp, integrator, p, ttmp)
    f(k₂, tmp, p, ttmp)
    @.. broadcast=false thread=thread tmp=uprev + halfdt * k₂
    stage_limiter!(tmp, integrator, p, ttmp)
    f(k₃, tmp, p, ttmp)
    @.. broadcast=false thread=thread tmp=uprev + dt * k₃
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k₄, tmp, p, t + dt)
    @.. broadcast=false thread=thread u=uprev + (dt / 6) * (2 * (k₂ + k₃) + (k₁ + k₄))
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    if integrator.opts.adaptive
        # Shampine Solving ODEs and DDEs with Residual Control Estimate
        k₅ = k
        _p = k₂
        pprime = k₃ # Alias some cache arrays
        # one(t) so that types are correct but unitless
        σ₁ = one(t) * (1 // 2) - sqrt(one(t) * 3) / 6
        σ₂ = one(t) * (1 // 2) + sqrt(one(t) * 3) / 6
        @.. broadcast=false thread=thread tmp=(1 - σ₁) * uprev + σ₁ * u +
                                              σ₁ * (σ₁ - 1) *
                                              ((1 - 2σ₁) * (u - uprev) +
                                               (σ₁ - 1) * dt * k₁ +
                                               σ₁ * dt * k₅)
        @.. broadcast=false thread=thread pprime=k₁ +
                                                 σ₁ *
                                                 (-4 * dt * k₁ - 2 * dt * k₅ - 6 * uprev +
                                                  σ₁ *
                                                  (3 * dt * k₁ + 3 * dt * k₅ + 6 * uprev -
                                                   6 * u) +
                                                  6 * u) / dt
        f(_p, tmp, p, t + σ₁ * dt)
        @.. broadcast=false thread=thread tmp=dt * (_p - pprime)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        e1 = integrator.opts.internalnorm(atmp, t)
        @.. broadcast=false thread=thread tmp=(1 - σ₂) * uprev + σ₂ * u +
                                              σ₂ * (σ₂ - 1) *
                                              ((1 - 2σ₂) * (u - uprev) +
                                               (σ₂ - 1) * dt * k₁ +
                                               σ₂ * dt * k₅)
        @.. broadcast=false thread=thread pprime=k₁ +
                                                 σ₂ *
                                                 (-4 * dt * k₁ - 2 * dt * k₅ - 6 * uprev +
                                                  σ₂ *
                                                  (3 * dt * k₁ + 3 * dt * k₅ + 6 * uprev -
                                                   6 * u) +
                                                  6 * u) / dt
        f(_p, tmp, p, t + σ₂ * dt)
        @.. broadcast=false thread=thread tmp=dt * (_p - pprime)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        e2 = integrator.opts.internalnorm(atmp, t)
        integrator.EEst = convert(typeof(one(t)), 2.1342) * max(e1, e2)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
end