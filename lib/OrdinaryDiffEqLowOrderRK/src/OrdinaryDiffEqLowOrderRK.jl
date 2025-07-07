module OrdinaryDiffEqLowOrderRK

import OrdinaryDiffEqCore: alg_order, isfsal, beta2_default, beta1_default,
                           alg_stability_size,
                           ssp_coefficient, OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqExponentialAlgorithm,
                           explicit_rk_docstring, generic_solver_docstring,
                           trivial_limiter!,
                           OrdinaryDiffEqAdaptiveAlgorithm,
                           unwrap_alg, @unpack, initialize!, perform_step!,
                           calculate_residuals,
                           calculate_residuals!, _ode_addsteps!, @OnDemandTableauExtract,
                           constvalue,
                           OrdinaryDiffEqMutableCache, uses_uprev,
                           OrdinaryDiffEqConstantCache, @fold,
                           @cache, CompiledFloats, alg_cache, CompositeAlgorithm,
                           copyat_or_push!,
                           AutoAlgSwitch, _ode_interpolant, _ode_interpolant!, full_cache,
                           accept_step_controller, DerivativeOrderNotPossibleError,
                           du_cache, u_cache, get_fsalfirstlast
using SciMLBase
import MuladdMacro: @muladd
import FastBroadcast: @..
import LinearAlgebra: norm
import RecursiveArrayTools: recursivefill!, recursive_unitless_bottom_eltype
import Static: False
using DiffEqBase: @def, @tight_loop_macros
import OrdinaryDiffEqCore

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("low_order_rk_caches.jl")
include("fixed_timestep_perform_step.jl")

export RK4

end
