

""" log_threshold
    A constant indicating when `lgamma` rather than `gamma` should be used in computation.
    """
const log_threshold = 100.0

""" logΓ(z)

    Returns the (single precision) natural log of Γ(z) where typeof(z) may be any subtype of Real, or of type Complex{T} where T is not multiprecision. Depending on the magnitude of z, SpecialFunctions.gamma or SpecialFunctions.lgamma is used. The two functions give comparable answers (i.e., exp(lgamma(z)) ≈ gamma(z) and lgamma(z) ≈ log(gamma(z))) when 100.0 < |z| < 150.0. According to documentation, gamma is more accurate for smaller magnitudes, lgamma for large.

    """
function logΓ(z::T) where T<:NotBigComplex
    return abs(z) > log_threshold ? lgamma(z) : log(gamma(z))
end


""" ζ(z)

    Alias for SpecialFunctions implementation, `zeta`, of Riemann's ζ function. Note: zeta can take real, but not complex, multiprecision arguments.
    """
ζ = zeta

""" ξ(s)

    Implementation of the Riemann xi function, ξ, using Riemann's zeta, ζ, and the gamma function, Γ, as implemented in Julia's SpecialFunctions package. Because zeta and gamma can take real, but not complex, multiprecision arguments, the same restrictions apply to ξ.
    """
function ξ(s::T) where {T<:Union{NotBigComplex, Real}}
    epsilon = eps(promote_type(typeof(real(s)),typeof(imag(s)),Float64))
    if abs(s-1.0) ≤ epsilon
        return (s/2)*π^(-s/2)*bigΓ(s/2) # (s-1)*ζ(s)→1 as s→1
    elseif abs(s) ≤ epsilon
        return π^(-s/2)*(s-1)*ζ(s) # s/2*Γ(s/2)→1 as s→0
    else
        return π^(-s/2)*(s/2)*bigΓ(s/2)*(s-1)*ζ(s)
    end
end

""" xi(s)

    Alias for Riemann's xi function, ξ(s). Note: because of restrictions in Julia's SpecialFunctions package, xi can take real, but not complex, multiprecision arguments.
    """
xi =  ξ

""" ΞRL(z)
    
    Implementation of the Riemann-Landau Xi function, Ξ, as ξ(1/2 + z*i).  Note: because of restrictions in Julia's SpecialFunctions package, ΞRL can take real, but not complex, multiprecision arguments.
     """
function ΞRL(z::T) where {T<:Union{NotBigComplex, Real}}
    return  ξ(1/2 + z*im)
end

""" XiRL(z)

    Alias for the Riemann-Landau Xi function, Ξ(z). Note: because of restrictions in Julia's SpecialFunctions package, Xi can take real, but not complex, multiprecision arguments.
    """
XiRL = ΞRL
    
