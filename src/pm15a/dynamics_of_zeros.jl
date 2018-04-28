""" log_r₀(n::Int, s::Number)

    We have, from equation (36) of the reference and Riemann-Siegel:

    r₀ₙ(s) = (1/8)×s(s-1)/2×√(π⁻ˢ)×Γ(s/2)×n⁻ˢ
    Here we implement log(r₀ₙ(s))
    """
function log_r₀(n::Int, s::Number)
    c₁ = s≈0 ? 0.0 : log(s/2) + lgamma(s/2) # log( s/2×Γ(s/2) ) with singularity at s≈0.0 removed.
    c₂ = -log(8) + log(s-1) - s*(log(π)/2 + n) # log( 1/8×(s-1)×√(π⁻ˢ)×n⁻ˢ )
    return c₁ + c₂
end

"""
    We have, from equation (39) of the reference:
    
    rₜₙ(s) = exp(-t/4×αₙ²)∫exp(-√(t)vαₙ)×r₀ₙ(s+√(t)v + t/2×αₙ)/√(π)×exp(-v²)dv
    
    where
    r₀ₙ(s) = (1/8)×s(s-1)/2×√(π⁻ˢ)×Γ(s/2)×n⁻ˢ
    
    and the default value of αₙ is α(s)-log(n) from equation (44), pp 14 of the reference.
    
    Here we return a function which computes the log of the nth integrand at v, including the external factor, 
    """
function log_integrand_rₜₙ(t::Real, s::Number, n::Int; αₙ::Number = α(s) - log(n))
    sign(imag(s)) == sign(imag(s+αₙ)) || error("ℑ(s) and ℑ(s+αₙ) must have the same sign.")
    c₀ = -log(8) -t/4*αₙ^2 - log(π)/2 # constant: log( 1/8×exp(-t/4×αₙ²)/√(π) )
    c₁ = -√(t)*αₙ              # multiplies v in log( (exp(-√(t)vαₙ) )
    c₂ = s+t/2*αₙ              # constant additive in argument s+√(t)v + t/2×αₙ
    c₃ = √(t)                  # mulitplier of v in same argument
    return function (v::Real) c₀ + c₁*v + log_r₀(n, c₂ + c₃*v) - v^2 end
end
