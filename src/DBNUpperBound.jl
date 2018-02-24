module DBNUpperBound

using QuadGK
using SpecialFunctions



# utility.jl and Ht.jl
export phi_decay, Φpm, Ht, ζ, zeta, Γ, gamma, ξ, xi, XiRL, ΞRL,H0
# KKL.jl
export ψKKL, ϕKKL, ΦKKL, ΞKKL
# Ht_asymptotics.jl
export A, B, C

include("types.jl")
include("special_fcts.jl") 
include("Ht_asymptotics.jl")
include("utility.jl")
include("Ht.jl")

# Submodules
include("KKL/KKL.jl")
using .KKL
include("numerics/Numerics.jl")
using .Numerics

end # module
