@doc explicit_rk_docstring("The canonical Runge-Kutta Order 4 method.
Uses a defect control for adaptive stepping using maximum error over the whole interval.",
    "RK4",
    references = "@article{shampine2005solving,
      title={Solving ODEs and DDEs with residual control},
      author={Shampine, LF},
      journal={Applied Numerical Mathematics},
      volume={52},
      number={1},
      pages={113--127},
      year={2005},
      publisher={Elsevier}
      }")
Base.@kwdef struct RK4{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RK4(stage_limiter!, step_limiter! = trivial_limiter!)
    RK4(stage_limiter!, step_limiter!, False())
end