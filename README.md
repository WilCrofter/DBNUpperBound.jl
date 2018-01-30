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

STATUS:

`phi_decay`, `Φ` (an alias for `phi_decay`), and `Ht`, are implemented and are compatible with Python versions according to unit tests. There's an expository Jupyter notebook in the /notebooks subdirectory.

TODO:

* Generalize `Φ` to take complex arguments with restricted imaginary parts.
* Implement $H_0$ in terms of ξ and ζ and compare with `Ht(0.0, ⋅)`.

In Julia 0.6, `Base.runtests()` seems to have bugs. Looking into this. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```