export ζ, ξ

ζ = zeta

""" ξ(s)

    Returns a multiprecision estimate of Riemann's xi, ξ, of the given argument. The implementation uses Riemann's zeta, ζ, and the log gamma, lgamma, as implemented in Julia's SpecialFunctions package. Because ζ and lgamma can take real, but not complex, multiprecision arguments, the same restrictions apply to ξ.

    Basic precision is determined by the ambient value of `precision(BigFloat)`, but calculations may result in expansion.
    """
function ξ(s::Number)
    s64= typeof(s) <: Complex ? Complex{Float64}(s) : s
    if s ≈ 1.0
        return (s/2)*big(π)^(-s/2)*big(e)^(lgamma(s64/2)) # (s-1)*ζ(s)→1 as s→1
    elseif s ≈ 0.0
        return big(π)^(-s/2)*(s-1)*ζ(s64) # s/2*Γ(s/2)→1 as s→0
    else
        return big(π)^(-s/2)*(s/2)*big(e)^(lgamma(s64/2))*(s-1)*ζ(s64)
    end
end
