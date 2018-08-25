# de Bruijn Newman Upper Bound

Julia code following [Polymath Project 15](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant). The original motivation was to explore the efficacy of Julia in this context.

Code covering sections 1 through 8 of the draft, [EFFECTIVE APPROXIMATION OF HEAT FLOW EVOLUTION OF THE RIEMANN XI FUNCTION, AND AN UPPER BOUND FOR THE DE BRUIJN-NEWMAN CONSTANT](https://github.com/teorth/dbn_upper_bound/blob/0a0f72933da63a9589c42e3d404f4fc5b88060a7/Writeup/debruijn.pdf) (PDF, commit 0a0f729), has been written, tested, and mostly documented. Current functionality is exported from a submodule, `DBNUpperBound.PM15a`. The package remains substantially behind [dbn_upper_bound](https://github.com/teorth/dbn_upper_bound) where the original computational work is being done.

Arbitrary precision arithmetic is used by default. (Julia wraps the GNU Multiple Precision Arithmetic Library (GMP) and the GNU MPFR Library. [Reference](https://docs.julialang.org/en/stable/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic-1)).

## Installation

This version requires Julia 1+. Earlier versions, which require Julia 0.6+, can of course be checked out from this repo as well. (Julia 1.0.0, being a major release, in not backward compatible with 0.x.x releases.) To install this version and its dependencies:

Note that the SpecialFunctions dependency (for gamma and Riemann's zeta) requires both gcc and gcc-fortran.

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

All unit tests pass:

```
$ julia test/runtests.jl
[ Info: Testing PM15a
[ Info: Testing introduction.jl
[ Info: Testing applying the fundamental solution to the heat equation.
[ Info: Testing elementary estimates
[ Info: Testing estimates for large x
[ Info: Testing new upper bound
[ Info: Testing PM15b
[ Info: Testing estimating sums
[ Info: DBNUpperBound passes tests.
```
