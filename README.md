# de Bruijn Newman Upper Bound

Julia code following, i.e., having made no original contributions to, [Polymath Problem 15](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).

Status: Rewriting and refactoring for better accord with the current draft, [EFFECTIVE APPROXIMATION OF HEAT FLOW EVOLUTION OF THE RIEMANN XI FUNCTION, AND AN UPPER BOUND FOR THE DE BRUIJN-NEWMAN CONSTANT](https://github.com/teorth/dbn_upper_bound/blob/master/Writeup/debruijn.pdf) (PDF). This means notebooks and docs are obsolete, and general chaos prevail. Things should improve shortly. 

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

In Julia 0.6, `Base.runtests()` seems to have bugs. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```