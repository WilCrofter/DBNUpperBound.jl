using DBNUpperBound
using Base.Test

function test_estimates_for_large_x()
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
        

    end # precision
end

test_estimates_for_large_x()
