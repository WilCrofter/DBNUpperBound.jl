# de Bruijn Newman Upper Bound

Julia code in support of [Polymath Problem 15](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant) recently proposed by Terence Tao. The official code repositories are [here](https://github.com/km-git-acc/dbn_upper_bound).

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

### NEWS:

Implemented and tested H<sub>t</sub> approximations except for "toys" which are to be done, and effective errors E1, E2, E3, and E3*. As with H<sub>t</sub>, the default quadrature routine is inadequate for E3.


### TODO:

* Document and test everything which lacks same.

### Issues:


In Julia 0.6, `Base.runtests()` seems to have bugs. Looking into this. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```