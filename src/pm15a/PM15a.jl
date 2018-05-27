""" PM15a
    
    A submodule based on a draft version of "EFFECTIVE APPROXIMATION OF HEAT FLOW EVOLUTION OF THE RIEMANN XI FUNCTION, AND AN UPPER BOUND FOR THE DEBRUIJN-NEWMAN CONSTANT" by D.H.J. POLYMATH.

    The draft version, referred to as "the reference" in PM15a docs, is that of commit 0a0f729 of the github.com/torth/dbn_upper_bound. The original may be found at url https://github.com/teorth/dbn_upper_bound/blob/0a0f72933da63a9589c42e3d404f4fc5b88060a7/Writeup/debruijn.pdf. There is also a copy in this module.
    """
module PM15a

using DBNUpperBound
using SpecialFunctions
include("utilities.jl")
include("introduction.jl")
include("notation.jl")
include("dynamics_of_zeros.jl")
include("applying_the_fundamental_solution_for_the_heat_equation.jl")
include("elementary_estimates.jl")
include("estimates_for_large_x.jl")
include("bounding_dirichlet_series.jl")
include("a_new_upper_bound.jl")


export δ₁, F


end
