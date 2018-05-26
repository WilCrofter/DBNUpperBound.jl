
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
            αsum = asum = lhs = rhs1 = rhs2 = 0.0
            for n in 1:N(t,x)
                αₙ = αterm(t,n,x,y)
                a₁ = αₙ[1]*αₙ[2]
                a₂ = Aterm(t,n,x,y) # = bᵗₙ(t,n)*big(e)^(-astar*log(n))
                # Test that the αterm and Aterm values are the same within tolerance.
                @test isapprox(a₁,a₂,rtol=rtol)
                # Accumulate their values to test against loss of precision thru summation
                αsum += a₁
                asum += a₂
                # Test that |bᵗₙ⋅n^(y-κ̄)/n^bstar| ≤ |bᵗₙ⋅n^y/n^bstar| + bᵗₙ⋅(n^|κ̄|-1)/n^(ℜ(bstar)-y)
                lh = abs(αₙ[1]/big(e)^(bstar*log(n)))
                rh1 = abs(bᵗₙ(t,n)*big(e)^((y-bstar)*log(n)))
                rh2 = bᵗₙ(t,n)*(big(e)^(abs(κ(t,x,y)))-1)/big(e)^(real(bstar)-y)
                @test lh ≤ rh1 + rh2
                # Accumulate for summation test
                lhs += lh
                rhs1 += rh1
                rhs2 += rh2
            end
            @test isapprox(αsum, asum, rtol=rtol)
            @test isapprox(γₜ(t,x,y)*αsum, γₜ(t,x,y)*asum, rtol=rtol)
            @test lhs ≤ rhs1 + rhs2
        end

        # Test inequality (78)
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            λ̂ = mollifiers(3,t)
            # Note that the right hand side may be negative
            @test abs(fₜ(t,x,y)[1]) ≥ bound78(t,x,y,λ̂)
        end


    end
    
end

test_bounding_dirichlet_series()
