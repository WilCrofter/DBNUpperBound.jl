module PM15a

using DBNUpperBound
using DBNUpperBound.Asymptotics

include("introduction.jl")
include("estimates_for_large_x.jl")


export M₀, logM₀, logM₀′, α, Mₜ, B₀
export region_5, in_region_5, bound20, bound21, bound22
export γ, κ, fₜ, ϵ̃, ϵₜₙ, s_star




end
