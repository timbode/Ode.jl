using Test
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "lib", "OrdinaryDiffEqLowOrderRK"))
Pkg.test("OrdinaryDiffEqLowOrderRK")
