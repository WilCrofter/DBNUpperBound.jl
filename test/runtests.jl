using DBNUpperBound
using QuadGK
using Base.MathConstants
using SpecialFunctions
using Test
using Random

@info("Testing PM15a")
include("test_introduction.jl")
include("test_applying_the_fundamental_solution.jl")
include("test_elementary_estimates.jl")
include("test_estimates_for_large_x.jl")
include("test_new_upper_bound.jl")
@info("Testing PM15b")
include("test_estimating_sums.jl")
@info("DBNUpperBound passes tests.")
