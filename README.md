# de Bruijn Newman Upper Bound

An attempt to provide Julia code in support of Terence Tao's [recently proposed Polymath Problem.](http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant) In particular, I'll try to maintain Julia code which roughly parallels [@KM's python repository](https://github.com/km-git-acc/dbn_upper_bound).

To install this package and its dependencies:

```
julia> Pkg.clone("https://github.com/WilCrofter/DBNUpperBound.jl.git")
```
To remain current:
```
julia> Pkg.update("DBNUpperBound")
```

### NEWS:

Implemented asymptotic approximations A<sup>eff</sup>, B<sup>eff</sup> as described [in the wiki](http://michaelnielsen.org/polymath1/index.php?title=Polymath15_test_problem) and lightly tested it against corresponding Python code.

Implemented asymptotic approximations A', B' as described [in the wiki](http://michaelnielsen.org/polymath1/index.php?title=Polymath15_test_problem) and lightly tested it against corresponding Python code.

Started using [Jupyter notebooks](http://jupyter.org/) when implementing formulae. Notebooks allow juxtaposition of code with mathematical notation copied and pasted (as images or LaTeX) from wiki or blog pages. In my case this seems to result in cleaner and less error-prone implementation.

Implemented asymptotic approximations [A, B, C](https://terrytao.wordpress.com/2018/02/12/polymath15-third-thread-computing-and-approximating-h_t/) as in Thread 3, and tested for expected behavior for x between 10<sup>5</sup> and 10<sup>15</sup>.

### TODO:

* Implement effective errors for A<sup>eff</sup> and B<sup>eff</sup>.
* Document and test everything which lacks same.
* Check Ht as implemented against asymptotic approximations based on A, B, C. More likely: reimplement using notebooks.

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