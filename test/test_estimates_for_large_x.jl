using DBNUpperBound
using Base.Test
include("Ht_ref_values.jl") # Htref, a Dictionary of Hₜ reference values provided by Anonymous

function test_estimates_for_large_x()
    info("Testing estimates for large x")
    setprecision(80)do
        # julia> UInt(Dates.value(Dates.now()))
        # 0x000039e6618e0dda
        srand(0x000039e6618e0dda)
        
        r5 = region_5(100)
        s5 = [s⁺(r5[i,2],r5[i,3]) for i in 1:100]
        
        # Inequality (42), pp 13
        # α′(s) = O≤(1/(2*ℑ(s)-6) whenever ℑ(s) > 3
        # Ensure ℑ(s) > 10 for subsequent tests.
        @test sum([imag(s)>10 for s in s5]) == length(s5)
        @test all([abs(α′(s)) ≤ 1/(2*imag(s)-6) for s in s5])

        # Inequality in 2nd line of proof, pp 14
        n5 = rand(1:10,length(s5))
        @test all([abs(r₀(n5[i],s5[i])) ≤ abs(M₀(s5[i])*n5[i]^(-s5[i]))*big(e)^(1/(6*abs(s5[i])-0.66)) for i in 1:length(s5)])

        # Identity, last line pp 14,
        # M₀(σ+iT)exp(t/4*αₙ²-(σ+iT)log(n)) = Mₜ(σ+iT)*bᵗₙ/n^(σ+iT+t/2*α(σ+iT))
        for i in 1:length(s5)
            s = s5[i]
            t = r5[i,1]
            n = n5[i]
            αₙ = α(s)-log(n)
            @test M₀(s)*big(e)^(t/4*αₙ^2-s*log(n)) ≈ Mₜ(t,s)*bᵗₙ(t,n)/n^(s+t/2*α(s))
        end

        # Inequality (45) pp 15, testing closed form and quadrature
        for i in 1:length(s5)
            t = r5[i,1]
            T = imag(s5[i])
            lhsq = quadgk((v)->big(e)^(t*v^2/(2*(T-3.08))-v^2)/√(pi),-10.0,10.0)
            lhse = √(1-big(t)/(2*(T-3.08)))
            rhs = big(e)^(t/(4(T-3.33)))
            @test lhsq[1] ≤ rhs+lhsq[2] && lhse ≤ rhs && lhse^2 ≤ rhs
        end

        # Test crude bound Hₜ(z) = A(z)+B(z)+O≤(EA(z)+EB(z)+EC₀(z))
        # The expression is taken to mean that
        #          |Hₜ(z)| ≤ |A(z)+B(z) + EA(z)+EB(z)+EC₀(z)|
        # which implies
        #          |Hₜ(z)| ≤ |A(z)+B(z)| + EA(z)+EB(z)+EC₀(z)
        # by the triangle inequality and the fact that the error terms are ≥0.
        #
        # Anonymous provided exact calculations of Hₜ for t=y=.4 and
        # x = 10, 30, 100, 300, 1000, 3000, 10000. Only values of x
        # exceeding 200 are in region 5.
        t=y=.4
        for x in [300, 1000, 3000, 10000]
            @test abs(Htref[x]) ≤ abs(A(t,x,y)+B(t,x,y)+EA(t,x,y)+EB(t,x,y)+EC₀(t,x,y))
        end

        # Test finer bound Hₜ(z) = A(z)+B(z)-C(z)+O≤(EA(z)+EB(z)+EC(z))
        for x in [300, 1000, 3000, 10000]
            @test abs(Htref[x]) ≤ abs(A(t,x,y)+B(t,x,y)-C(t,x,y)+EA(t,x,y)+EB(t,x,y)+EC(t,x,y))
        end

        # Test bounds for ϵₜₙ, ẽA, ẽB
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            n = n5[i]
            s = s⁺(x,y)
            T = imag(s)
            @test (t^2/8*abs(α(s)-log(n))^2+t/4+1/6)/(T-3.33) ≤ (t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66)
            T = imag(1-s') # same as imag(s) of course
            @test (t^2/8*abs(α(1-s')-log(n))^2+t/4+1/6)/(T-3.33) ≤ (t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66)
            @test abs(ϵₜₙ(t,n,s)) ≤ abs(ϵ̃ₜₙ(t,n,s))
            @test abs(ϵₜₙ(t,n,1-s')) ≤ abs(ϵ̃ₜₙ(t,n,1-s'))
            # @test eA(t,x,y) ≤ ẽA(t,x,y)
            @test abs(eB(t,x,y)) ≤ abs(ẽB(t,x,y))
        end

    end # precision
end

test_estimates_for_large_x()
