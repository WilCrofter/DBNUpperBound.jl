# de Bruijn Newman Upper Bound

An attempt to provide Julia code in support of Terence Tao's [recently proposed Polymath Problem.](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant) In particular, I'll try to maintain Julia code which roughly parallels [@KM's python repository](https://github.com/km-git-acc/dbn_upper_bound).

The main point of Julia in this case is, of course, performance.

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

### STATUS:

`phi_decay`, `Φ` (an alias for `phi_decay`), and `Ht`, are implemented and are compatible with Python versions according to unit tests. `Φ` allows complex arguments with imaginary parts less than π/8 in absolute value. There's an expository Jupyter notebook in the /notebooks subdirectory.

Re-exported SpecialFunctions zeta, gamma, and exported respective aliases, ζ, and Γ. As implemented by the SpecialFunctions package, multiprecision complex arguments are disallowed. Multiprecision reals are OK.

Implemented ξ, its alias xi, Ξ, its alias Xi, and H0. Multiprecision complex arguments to ξ, xi, Ξ, Xi are disallowed due to dependence on SpecialFunctions zeta and gamma. Multiprecision arguments to H0 are disallowed.

### TODO:

* Unit tests for ξ, Ξ, and for H0 against Ht(0.0, ⋅).
* Expand docs to include optional parameters.

### Issues:

In Julia 0.6, `Base.runtests()` seems to have bugs. Looking into this. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```