""" PM15b
    
    A submodule based primarily on sections 7 & 8 of the July 8, 2018 (2018-7-8) draft of "EFFECTIVE APPROXIMATION OF HEAT FLOW EVOLUTION OF THE RIEMANN XI FUNCTION, AND AN UPPER BOUND FOR THE DEBRUIJN-NEWMAN CONSTANT" by D.H.J. POLYMATH.

    The 2018-7-8 draft, referred to as "the reference" in PM15b docs, is that of commit 3f2f5bd of the github.com/torth/dbn_upper_bound. The original may be found at url https://github.com/teorth/dbn_upper_bound/blob/3f2f5bd6949df68b9ba8a77bc3020c3453548f97/Writeup/debruijn.pdf. There is also a copy in this module. (The PDF uses git Large File Storage, though that should be transparent to forks and clones.)
    """
module PM15b

using DBNUpperBound
using SpecialFunctions

include("estimating_sums.jl")


end
