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

export Hₜ, ζ
export s⁺, M₀, logM₀, logM₀′, α, Mₜ, B₀
export region_5, in_region_5, bound20, bound21, bound22, bound23, bound24
export γₜ, κ, fₜ, s_star

export complex_power, gaussian_identity # temporarily, anyway
export r₀, rₜₙ_integrand, rₜₙ_by_integration

# estimates for large x
export ϵₜₙ, rₜₙ, RtN, A, B, C, C₀, EA, EB, EC
export ϵ̃, eA, eB, eC, eC0

export δ₁, F



    
end # module
