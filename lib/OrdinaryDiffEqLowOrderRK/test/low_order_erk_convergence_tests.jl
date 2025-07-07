# This definitely needs cleaning
using OrdinaryDiffEqLowOrderRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random
using DiffEqBase
using SciMLBase

Random.seed!(100)

## Convergence Testing
dts1 = 1 .// 2 .^ (9:-1:5)
dts2 = 1 .// 2 .^ (7:-1:3)
dts3 = 1 .// 2 .^ (12:-1:7)
dts4 = 1 .// 2 .^ (5:-1:3)
dts5 = 1 .// 2 .^ (3:-1:1)
dts6 = 1 .// 10 .^ (5:-1:1)
testTol = 0.2

f = (u, p, t) -> sin(u)
prob_ode_nonlinear = ODEProblem(
    ODEFunction(f;
        analytic = (u0, p, t) -> 2 * acot(exp(-t) *
                                          cot(0.5))), 1.0,
    (0.0, 0.5))

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear)[i]
    dts = 1 .// 2 .^ (8:-1:4)
    @info "Very low order"
    sim3 = test_convergence(dts, prob, RK4())
    @test sim3.𝒪est[:l∞]≈4 atol=testTol
   end

# Re-export the most common user-facing types/functions
export ODEProblem, ODEFunction, solve
