""" Ht_term(t, u, n)

    Returns the value of the nth term of the infinite series, Φ(u)*exp(t*u^2), where Φ(u) is  defined as on the PolyMath 15 page http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant. This is written separately from the corresponding term (Φpm(u)) of Φ(u) itself to avoid numerical problems with large positive and negative exponents.
    """
function Ht_term(t::T1, u::T2, n::Int) where {T1<:Real, T2<:Real}
    x = t*u^2-π*n^2*exp(4*u)
    return (2*π^2*n^4*exp(9*u+x)-3*π*n^2*exp(5*u+x))
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

""" ζ(z)

    Alias for SpecialFunctions implementation, `zeta`, of Riemann's ζ function. Note: zeta allows BigFloat and BigInt arguments, but not Complex{BigFloat} or Complex{BigInt}.
    """
ζ = zeta

""" Γ(z)

    Alias for SpecialFunctions implemention, `gamma` of the gamma function.  Note: gamma allows BigFloat and BigInt arguments, but not Complex{BigFloat} or Complex{BigInt}.
    """
Γ = gamma


""" NotBigReal

    A convenient type for trapping certain arguments. SpecialFunctions `gamma` and `zeta` can take real, but not complex, BigFloats.
    """
const NotBigReal = Union{Signed, Rational, Float64, Float32, Float16}

""" NotBigComplex

    A convenient type for trapping certain arguments. SpecialFunctions `gamma` and `zeta` can take real, but not complex, BigFloats.
    """
const NotBigComplex = Complex{T} where {T <: NotBigReal}

""" ξ(s)

    Implementation of the Riemann xi function, ξ, using Riemann's zeta, ζ, and the gamma function, Γ, as implemented in Julia's SpecialFunctions package. Because zeta and gamma can take real, but not complex, multiprecision arguments, the same restrictions apply to ξ.
    """
function ξ(s::T) where {T<:Union{NotBigComplex, Real}}
    epsilon = eps(promote_type(typeof(real(s)),typeof(imag(s)),Float64))
    if abs(s-1.0) ≤ epsilon
        return (s/2)*π^(-s/2)*Γ(s/2) # (s-1)*ζ(s)→1 as s→1
    elseif abs(s) ≤ epsilon
        return π^(-s/2)*(s-1)*ζ(s) # s/2*Γ(s/2)→1 as s→0
    else
        return π^(-s/2)*(s/2)*Γ(s/2)*(s-1)*ζ(s)
    end
end

""" xi(s)

    Alias for Riemann's xi function, ξ(s). Note: because of restrictions in Julia's SpecialFunctions package, xi can take real, but not complex, multiprecision arguments.
    """
xi =  ξ

""" Ξ(z)
    
    Implementation of the Riemann-Landau Xi function, Ξ, as ξ(1/2 + z*i).  Note: because of restrictions in Julia's SpecialFunctions package, Ξ can take real, but not complex, multiprecision arguments.
     """
function Ξ(z::T) where {T<:Union{NotBigComplex, Real}}
    return  ξ(1/2 + z*im)
end

""" Xi(z)

    Alias for the Riemann-Landau Xi function, Ξ(z). Note: because of restrictions in Julia's SpecialFunctions package, Xi can take real, but not complex, multiprecision arguments.
    """
Xi = Ξ
    
""" H0(z)

    Implementation of H0(z) as (1/8)*ξ(1/2+z*im/2). See https://en.wikipedia.org/wiki/Riemann_Xi_function. Compare with Ht(0.0,z). Note: because of restrictions in Julia's SpecialFunctions package, H0 cannot take multiprecision arguments.
    """
function H0(z::T) where {T<:Union{NotBigReal,NotBigComplex}}
    return  ξ((1+z*im)/2)/8
end
    
