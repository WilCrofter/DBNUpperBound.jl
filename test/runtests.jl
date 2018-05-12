using DBNUpperBound
using QuadGK
using SpecialFunctions
using Base.Test

info("Testing PM15a")
include("test_introduction.jl")
include("test_elementary_estimates.jl")
include("test_estimates_for_large_x.jl")
info("PM15a passes tests.")
