module DBNUpperBound

using QuadGK
using SpecialFunctions

Γ = SpecialFunctions.gamma

export bigexp, bigify

# utility.jl and Ht.jl
# export phi_decay, Φpm, Ht, ζ, zeta, Γ, gamma, ξ, xi, XiRL, ΞRL,H0
# KKL.jl
# export ψKKL, ϕKKL, ΦKKL, ΞKKL
# Ht_asymptotics.jl
# export A, B, C

# include("types.jl")
# include("special_fcts.jl") 
# include("utility.jl")
# include("Ht.jl")

# Submodules
# include("KKL/KKL.jl")
# using .KKL
# include("asymptotics/asymptotics.jl")
# using .Asymptotics
include("pm15a/PM15a.jl")
using .PM15a

export s⁺, M₀, logM₀, logM₀′, α, Mₜ, B₀
export region_5, in_region_5, bound20, bound21, bound22, bound23, bound24
export γₜ, κ, fₜ, ϵ̃, ϵₜₙ, s_star

export complex_power, gaussian_identity # temporarily, anyway
export log_r₀, log_integrand_rₜₙ

function bigify(s::T) where {T<:Number}
    return imag(s)==0 ? big(real(s)) : big(real(s))+im*big(imag(s))
end

function bigexp(s::T) where {T<:Number}
    ans = exp(real(s))
    ans = exp(im*imag(s))*(isfinite(ans) && !(ans == 0) ? ans : exp(bigify(real(s))))
    return imag(ans) ≈ 0.0 ? real(ans) : ans
end

end # module
