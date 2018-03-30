# de Bruijn Newman Upper Bound

Julia code for [Polymath Problem 15](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant). This repository generally lags the [official repository](https://github.com/km-git-acc/dbn_upper_bound)--my interests are somewhat tangential to the mainstream and analytic number theory is not my field.

The notebook directory above is probably the best indicator of what's going on here.

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

### TODO:

* Document and test everything which lacks same.

### Issues:


In Julia 0.6, `Base.runtests()` seems to have bugs. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```