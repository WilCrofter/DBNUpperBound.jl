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

`phi_decay`, `Φpm` (an alias for `phi_decay`), and `Ht`, are implemented and are compatible with Python versions according to unit tests. `Φpm` allows complex arguments with imaginary parts less than π/8 in absolute value. There's an expository Jupyter notebook in the /notebooks subdirectory.

Re-exported SpecialFunctions zeta, gamma, and exported respective aliases, ζ, and Γ. As implemented by the SpecialFunctions package, multiprecision complex arguments are disallowed. Multiprecision reals are OK.

Implemented ξ, its alias xi, Riemann & Landau's Ξ, ΞRL, its alias XiRL, and H0. Multiprecision complex arguments to ξ, xi, ΞRL, XiRL are disallowed due to dependence on SpecialFunctions zeta and gamma. Multiprecision arguments to H0 are disallowed.

Implemented ψ, ϕ, Φ, Ξ<sub>λ</sub> as described in [KKL2009](https://www.sciencedirect.com/science/article/pii/S0001870809001133).  

Tested approximate equalities Ht(0,t)≈1/2*Ξ(0.0,t/2) and Φpm(t)≈1/2*ΦKKL(2*t) which can be deduced from definitions.

### TODO:

* Implement K_tθ and I_tθ, -π/8 < θ < π/8, as described [in first thread](https://terrytao.wordpress.com/2018/01/27/polymath15-first-thread-computing-h_t-asymptotics-and-dynamics-of-zeroes/).
* Expand docs to include optional parameters.
* Document ψ, ϕ, ΦKKL, Ξ(λ,z).

### Issues:

`Φpm(z)` as computed is not an even function within numerical tolerance. Looking into it.

In Julia 0.6, `Base.runtests()` seems to have bugs. Looking into this. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```