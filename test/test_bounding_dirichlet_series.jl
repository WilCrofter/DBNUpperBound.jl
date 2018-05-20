
using DBNUpperBound
using Base.Test
if !isdefined(:Href)
    include("Ht_ref_values.jl") # Htref, a Dictionary of Hₜ reference values provided by Anonymous
end

function test_bounding_dirichlet_series()
    info("Testing pm15a/bounding_dirichlet_series.jl")

    setprecision(127)do
        
        # test mollifiers
        λ = mollifiers(3,.2)
        @test λ[1] == 1
        @test λ[2] == -bᵗₙ(.2,2)
        @test λ[3] == -bᵗₙ(.2,3)
        @test λ[6] == bᵗₙ(.2,2)*bᵗₙ(.2,3)

        # julia> hash(now())
        # 0xe1975fbaad83093b
        srand(0xe1975fbaad83093b)
        r5 = region_5(10)

        #=
        We have that fₜ = A/B₀ + B/B₀ by (14) pp. 4, the definitions of A and B in Cor 6.4 pp. 22, 
        the definition of B₀ in (11) pp. 3, and the definition of γ in (16) pp. 4.
        Moreover fₜ = A/B₀ + B/B₀ describes how it fₜ is computed in introduction.jl. 

        Computationally, then, we should have |fₜ| ≥ ||B/B₀|-|A/B₀||.
        =#
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            @test abs(fₜ(t,x,y)[1]) ≥ abs( abs(B(t,x,y)[1]/B₀(t,x,y)) - abs(A(t,x,y)[1]/B₀(t,x,y)) )
        end
        
        
        
    end
    
end

test_bounding_dirichlet_series()
