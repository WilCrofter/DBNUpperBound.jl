export r₀, rₜₙ_integrand, rₜₙ_by_integration

""" r₀(n::Int, s::Number)

    We have, from equation (36) of the reference and Riemann-Siegel:

    r₀ₙ(s) = (1/8)×s(s-1)/2×√(π⁻ˢ)×Γ(s/2)×n⁻ˢ
           = (1/8)×(s-1)×√(π⁻ˢ)×n⁻ˢ×s/2Γ(s/2)
    """
function r₀(n::Int, s::Number)
    c₁ = s≈0 ? 1.0 : s/2*big(e)^(lgamma(s/2)) # s/2×Γ(s/2) with singularity at s≈0.0 removed.
    c₂ = 1/8*(s-1)*big(e)^(-s*(log(π)/2 + n)) # 1/8×(s-1)×√(π⁻ˢ)×n⁻ˢ
    return c₁*c₂
end

"""
    We have, from equation (39) of the reference:
    
    rₜₙ(s) = exp(-t/4×αₙ²)∫exp(-√(t)vαₙ)×r₀ₙ(s+√(t)v + t/2×αₙ)/√(π)×exp(-v²)dv
    
    where
    r₀ₙ(s) = (1/8)×s(s-1)/2×√(π⁻ˢ)×Γ(s/2)×n⁻ˢ
    
    and the default value of αₙ is α(s)-log(n) from equation (44), pp 14 of the reference.
    
    Here we return a function which computes the nth integrand at v, including the external factor, 
    """
function rₜₙ_integrand(t::Real, s::Number, n::Int; αₙ::Number = α(s) - log(n))
    sign(imag(s)) == sign(imag(s+αₙ)) || error("ℑ(s) and ℑ(s+αₙ) must have the same sign.")
    c₀ = big(e)^(-t/4*αₙ^2)     # external factor
    c₁ = -√(t)*αₙ              # multiplies v in exp(-√(t)vαₙ)
    c₂ = s+t/2*αₙ              # constant additive in argument s+√(t)v + t/2×αₙ
    c₃ = √(t)                  # mulitplier of v in same argument
    return function (v::Real) c₀*big(e)^(c₁*v-v^2)*r₀(n, c₂ + c₃*v) end
end



function rₜₙ_by_integration(t::Real, s::Number, n::Int; αₙ::Number = α(s) - log(n),
                            lower_limit::Float64=-10.0, upper_limit::Float64=10.0)
    g = rₜₙ_integrand(t,s,n,αₙ=αₙ)
    return quadgk(g,upper_limit,0.0,lower_limit)
end
