""" Ht_term(t, u, n)

    Returns the value of the nth term of the infinite series, Φ(u)*exp(t*u^2), where Φ(u) is  defined as on the PolyMath 15 page http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant. This is written separately from the corresponding term (Φpm(u)) of Φ(u) itself to avoid numerical problems with large positive and negative exponents.
    """
function Ht_term(t::T1, u::T2, n::Int) where {T1<:Real, T2<:Real}
    x = t*u^2-π*n^2*exp(4*u)
    return (2*π^2*n^4*exp(9*u+x)-3*π*n^2*exp(5*u+x))
end

"""
    If n is real, it is rounded down to the nearest integer.
    """
function Ht_term(t::T1, u::T2, v::Float64) where {T1<:Real, T2<:Real}
    Ht_term(t,u,floor(Int,v))
end

"""
    Sum Ht terms to nmax and estimate the tail and tail error by quadrature. Tail and tail errors are generally VERY small.
    """
function Ht_sum(t, u; nmax=10)
    s = 0.0
    for n in 1:nmax
        s += Ht_term(t, u, n)
    end
    tail, err = quadgk((v::Float64)->Ht_term(t,u,v),nmax+1.5,Inf,abstol=eps(Float64))
    s+tail, err
end

""" Ht_integrand(t, u, z; nmax=100)
    
    Returns a partial summation of the infinite series, Φ(u)*exp(t*u^2), to n_max (default 100) terms, where Φ(u) is  defined as on the PolyMath 15 page http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant.
    """
function Ht_integrand(t::T1, u::T2, z::T3; n_max=100) where {T1<:Real, T2<:Real, T3<:Number}
    ans = convert(promote_type(Complex{T1},Complex{T2},T3),0)
    for n in 1:n_max
        ans += Ht_term(t,u,n)
    end
    return ans*cos(z*u)
end

""" Ht(t, z; n_max=100)
    
    Returns a 2-tuple of an approximate value of Ht(z) and an estimated error of approximation. Ht(z) is the integral from 0 to ∞ of Φ(u)*exp(t*u^2)*cos(uz), where the first two factors are approximated by truncated series of n_max (default 100) terms. (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).
    """
function Ht(t::T1, z::T2; n_max::Int=100, upper_limit::T1=10.0, abstol=eps(T1), maxevals=10^7) where {T1<:Real,T2<:Number}
    return quadgk((u)-> Ht_integrand(t,u,z; n_max=n_max), 0.0, upper_limit, abstol=abstol, maxevals=maxevals)
end

""" H0(z)

    Implementation of H0(z) as (1/8)*ξ(1/2+z*im/2). See https://en.wikipedia.org/wiki/Riemann_Xi_function. Compare with Ht(0.0,z). Note: because of restrictions in Julia's SpecialFunctions package, H0 cannot take multiprecision arguments.
    """
function H0(z::T) where {T<:Union{NotBigReal,NotBigComplex}}
    return  ξ((1+z*im)/2)/8
end
    
