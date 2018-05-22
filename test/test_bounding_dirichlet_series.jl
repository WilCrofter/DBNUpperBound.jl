
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
        r5 = region_5(25)

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
        
        # Test inequality (77)
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            @test abs(fₜ(t,x,y)[1]) ≥ bound77(t,x,y)
        end

        # Test numerical consistency of Aterm and αterm
        rtol = (precision(BigFloat)-1)/log(2,10) # essentially eps(BigFloat)
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            s = s⁺(x,y)
            astar = s+t/2*α(s)
            bstar = 1-s+t/2*α(1-s)
            # is bstar-astar-(y-κ̄) imaginary within tolerance?
            @test isapprox(abs(big(e)^(bstar-astar-(y-κ(t,x,y)'))),1,rtol=rtol)
            # is bstar-astar ≈ (y-κ̄) + ix+t⋅(α(s̄)-α(s))/2 within tolerance?
            @test isapprox(bstar-astar, y-κ(t,x,y)'+im*x+t/2*(α(s')-α(s)),rtol=rtol)
            # αterm returns two numbers whose product should be bᵗₙ/n^astar within tolerance
            αsum = asum = 0.0
            for n in 1:N(t,x)
                αₙ = αterm(t,n,x,y)
                a₁ = αₙ[1]*αₙ[2]
                a₂ = Aterm(t,n,x,y) # = bᵗₙ(t,n)*big(e)^(-astar*log(n))
                @test isapprox(a₁,a₂ , rtol=rtol)
                αsum += a₁
                asum += a₂
            end
            @test isapprox(αsum, asum, rtol=rtol)
            @test isapprox(γₜ(t,x,y)*αsum, γₜ(t,x,y)*asum, rtol=rtol)
        end

    end
    
end

test_bounding_dirichlet_series()
