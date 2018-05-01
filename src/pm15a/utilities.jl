function bigify(s::T) where {T<:Number}
    return imag(s)==0 ? big(real(s)) : big(real(s))+im*big(imag(s))
end

function bigexp(s::T) where {T<:Number}
    ans = exp(real(s))
    ans = exp(im*imag(s))*(isfinite(ans) && !(ans == 0) ? ans : exp(bigify(real(s))))
    return imag(ans) ≈ 0.0 ? real(ans) : ans
end

function Hₜ(t::Real, x::Real, y::Real; lower_limit::Real=-10.0, upper_limit::Real=+10.0) 
    s=s⁺(x,y)
    srt = √(t)
    q(v) = ξ(s+srt*v)*bigexp(-v^2)/8/√(π)
    return quadgk(q,lower_limit,0.0,upper_limit)
end

ζ = zeta

""" ξ(s)

    Implementation of the Riemann xi function, ξ, using Riemann's zeta, ζ, and the log gamma, lgamma, as implemented in Julia's SpecialFunctions package. Because zeta and gamma can take real, but not complex, multiprecision arguments, the same restrictions apply to ξ.
    """
function ξ(s::T) where T<:Union{Float64, Complex{Float64}}
    epsilon = eps(promote_type(typeof(real(s)),typeof(imag(s)),Float64))
    if abs(s-1.0) ≤ epsilon
        return (s/2)*π^(-s/2)*bigexp(lgamma(s/2)) # (s-1)*ζ(s)→1 as s→1
    elseif abs(s) ≤ epsilon
        return π^(-s/2)*(s-1)*ζ(s) # s/2*Γ(s/2)→1 as s→0
    else
        return π^(-s/2)*(s/2)*bigexp(lgamma(s/2))*(s-1)*ζ(s)
    end
end
