# de Bruijn Newman Upper Bound

Julia code following, i.e., having made no original contributions to, [Polymath Problem 15](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).

Status: Refactored for clarity, better numerics, and better accord with the current draft, [EFFECTIVE APPROXIMATION OF HEAT FLOW EVOLUTION OF THE RIEMANN XI FUNCTION, AND AN UPPER BOUND FOR THE DE BRUIJN-NEWMAN CONSTANT](https://github.com/teorth/dbn_upper_bound/blob/master/Writeup/debruijn.pdf) (PDF). Still substantially behind [dbn_upper_bound](https://github.com/teorth/dbn_upper_bound) and still a work in progress. There is very little documentation at the moment.

Default use of arbitrary precision arithmetic is a major difference w.r.t. earlier code. (Julia wraps the GNU Multiple Precision Arithmetic Library (GMP) and the GNU MPFR Library. [Reference](https://docs.julialang.org/en/stable/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic-1)). 

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

In Julia 0.6, `Base.runtests()` seems to fail under certain circumstances. If so, the following alternative should work:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```