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
include("bounding_dirichlet_series.jl")
include("a_new_upper_bound.jl")


export δ₁, F


end
