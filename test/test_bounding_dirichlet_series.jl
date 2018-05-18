
using DBNUpperBound
using Base.Test

function test_bounding_dirichlet_series()
    info("Testing pm15a/bounding_dirichlet_series.jl")

    setprecision(80)do
        # test mollifiers
        λ = mollifiers(3,.2)
        @test λ[1] == 1
        @test λ[2] == -bᵗₙ(.2,2)
        @test λ[3] == -bᵗₙ(.2,3)
        @test λ[6] == bᵗₙ(.2,2)*bᵗₙ(.2,3)
    end
end

test_bounding_dirichlet_series()
