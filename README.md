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

Implemented asymptotic approximations [A, B, C](https://terrytao.wordpress.com/2018/02/12/polymath15-third-thread-computing-and-approximating-h_t/) as in Thread 3, and tested for expected behavior for x between 10<sup>5</sup> and 10<sup>15</sup>.

### TODO:

* Check Ht as implemented against asymptotic approximations based on A, B, C.
* Implement logA, logB, logC in attempt to increase range?
* Expand docs to include optional parameters.
* Document ψ, ϕ, ΦKKL, Ξ(λ,z).

### Issues:


As implemented here, asymptotic approximations, A+B-C, disagree with Ht by orders of magnitude for large values of the real part of z. In some sense this is not surprising since for large real(z) Ht's reported quadrature error is larger than its estimate of the integral, but it's suspicious enough to warrant study. Ht, as implemented here, does agree with its Python counterpart for smaller values of real(z). (See `test/test_Ht.jl`.) Example:
```
julia> D.Ht(0.4,1e3+.4im)
(-1.1686937089741209e-18 + 2.5500530718937875e-20im, 2.2322557270968237e-16)

julia> z = 1e3 + .4im
1000.0 + 0.4im

julia> s=(1+im*z)/2
0.3 + 500.0im

julia> M=N=floor(Int,sqrt(imag(s)/(2*π)))
8

julia> A(0.4,N,s)+B(.4,M,s)-C(.4,N,M,s)
4.345878346665737394544694943915155845207144437164468917444559157180353968548644e-167 + 2.820807643087679411425141633331924224210463568073841825703857716030970761075199e-167im
```

`Φpm(z)` as computed is not an even function within numerical tolerance. Looking into it.


In Julia 0.6, `Base.runtests()` seems to have bugs. Looking into this. Meanwhile, the following alternative may serve:

```
julia> cd(Pkg.dir("DBNUpperBound"))
julia> include("test/runtests.jl")
Test Passed
```