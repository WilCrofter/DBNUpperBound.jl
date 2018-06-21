module DBNUpperBound

using QuadGK
using SpecialFunctions
using Primes

include("pm15a/PM15a.jl")
using .PM15a

# pm15a/utilities.jl
export ζ, ξ

# pm15a/applying_the_fundamental_solution_for_the_heat_equation.jl
export Hₜ

# pm15a/dynamics_of_zeros.jl
export r₀, rₜₙ_integrand, rₜₙ_by_integration

# pm15a/introduction.jl
export s⁺, M₀, logM₀, logM₀′, α, α′, Mₜ, B₀, bᵗₙ, logbᵗₙ
export region_5, in_region_5
export γₜ, κ, fₜ
export N, H̃, H̃₂

# pm15a/notation
export complex_power, gaussian_identity # temporarily, anyway

# pm15a/estimates_for_large_x
export ϵₜₙ, ϵ̃ₜₙ, ϵ̃, rₜₙ, RtN, A, B, C, EA, EB, EC, EC₀
export eA, eB, eC, eC₀
export Aterm, Bterm

# pm15a/bounding_dirichlet_series
export Dirichlet_convolution, mollifiers, αterm
export bound77, bound78

# pm15a/a_new_upper_bound.jl
export δ₁, δ₁_ub, bound23, bound23a, util_23b, bound23b, bound85, boundsδfₜ

# pm15a/estimating_sums.jl
export NᵢRange, FN, FN_tail, EBound

    
end # module
