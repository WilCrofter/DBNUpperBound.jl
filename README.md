# de Bruijn Newman Upper Bound

Julia code following [Polymath Project 15](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant). The original motivation was to explore the efficacy of Julia in this context.

Code covering sections 1 through 8 of the draft, [EFFECTIVE APPROXIMATION OF HEAT FLOW EVOLUTION OF THE RIEMANN XI FUNCTION, AND AN UPPER BOUND FOR THE DE BRUIJN-NEWMAN CONSTANT](https://github.com/teorth/dbn_upper_bound/blob/0a0f72933da63a9589c42e3d404f4fc5b88060a7/Writeup/debruijn.pdf) (PDF, commit 0a0f729), has been written, tested, and mostly documented. Current functionality is exported from a submodule, `DBNUpperBound.PM15a`. The package remains substantially behind [dbn_upper_bound](https://github.com/teorth/dbn_upper_bound) where the original computational work is being done.

Default use of arbitrary precision arithmetic is a major difference w.r.t. earlier code. (Julia wraps the GNU Multiple Precision Arithmetic Library (GMP) and the GNU MPFR Library. [Reference](https://docs.julialang.org/en/stable/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic-1)).

## Installation

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

In Julia 0.6+, `Base.runtests()` seems to fail under certain circumstances. If so, the following alternative should work:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
INFO: Testing PM15a
INFO: Testing introduction.jl
INFO: Testing applying the fundamental solution to the heat equation.
INFO: Testing elementary estimates
INFO: Testing estimates for large x
INFO: Testing new upper bound
INFO: PM15a passes tests.
```