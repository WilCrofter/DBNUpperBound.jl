
using DBNUpperBound
using QuadGK
using Base.MathConstants
using SpecialFunctions
using Test
using Random

# passes but takes time #

function test_new_upper_bound()
    @info("Testing new upper bound")
    
    setprecision(80)do
        # julia> hash(now())
        # 0xba15e93cb10b7401
        Random.seed!(0xba15e93cb10b7401)

        r5 = region_5(25)

        # test that eA+eB ≤ bound23
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            @test eA(t,x,y)+eB(t,x,y) ≤ bound23(t,x,y)
        end

        # test that δ₁ satisfies its upper bound, that bound23 ≤ bound23a,
        # that  |γ||N|^(|κ|+y)⋅(e^δ₁(t₀,x₀)-1) ≤ util_23b,
        # and that eC₀ ≤ bound85
        t₀ = maximum(r5[:,1])
        x₀ = minimum(r5[:,2])
        y₀ = minimum(r5[:,3])
        N₀ = N(t₀,x₀)
        δ_ub = δ₁(t₀,x₀)
        u23b = util_23b(t₀, x₀, 1.0)
        for i in 1:size(r5,1)
            t,x,y = r5[i,1],r5[i,2],r5[i,3]
            @test δ₁(t,x) ≤ δ_ub
            @test bound23(t,x,y) ≤ bound23a(t₀,x₀,t,x,y)
            @test abs(γₜ(t,x,y)) * big(e)^((abs(κ(t,x,y))+y)*log(N(t,x))) *
                (big(e)^δ_ub-1) ≤ u23b
            @test bound23a(t₀,x₀,t,x,y) ≤ bound23b(t₀,x₀,t,x,y)
            @test eC₀(t,x,y) ≤ bound85(t, x₀, y₀, N₀)
        end

    end
end

test_new_upper_bound()
