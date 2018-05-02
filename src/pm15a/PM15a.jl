module PM15a

using DBNUpperBound
# using DBNUpperBound.Asymptotics
using SpecialFunctions
include("utilities.jl")
include("introduction.jl")
include("notation.jl")
include("dynamics_of_zeros.jl")
include("applying_the_fundamental_solution_for_the_heat_equation.jl")
include("elementary_estimates.jl")
include("estimates_for_large_x.jl")


export s⁺, M₀, logM₀, logM₀′, α, Mₜ, B₀
export region_5, in_region_5, bound20, bound21, bound22, bound23, bound24
export γₜ, κ, fₜ, ϵ̃, ϵₜₙ, s_star

export complex_power, gaussian_identity # temporarily, anyway
export log_r₀, log_integrand_rₜₙ

export Hₜ, ζ




end
