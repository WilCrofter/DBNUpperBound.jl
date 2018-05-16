module DBNUpperBound

using QuadGK
using SpecialFunctions

include("pm15a/PM15a.jl")
using .PM15a

# pm15a/utilities.jl
export ζ, ξ

# pm15a/applying_the_fundamental_solution_for_the_heat_equation.jl
export Hₜ

# pm15a/introduction.jl
export s⁺, M₀, logM₀, logM₀′, α, α′, Mₜ, B₀, bᵗₙ, logbᵗₙ
export region_5, in_region_5
export γₜ, κ, fₜ
export N, H̃, H̃₂


# pm15a/notation
export complex_power, gaussian_identity # temporarily, anyway

# pm15a/applying_the_fundamental ...
export r₀, rₜₙ_integrand, rₜₙ_by_integration

# pm15a/estimates_for_large_x
export ϵₜₙ, ϵ̃ₜₙ, ϵ̃, rₜₙ, RtN, A, B, C
#export C₀, EA, EB, EC, EC₀
#export ϵ̃, eA, eB, eC, eC0
#export ϵ̃ₜₙ, ẽA, ẽB

export δ₁, F



    
end # module
