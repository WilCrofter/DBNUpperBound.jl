
using DBNUpperBound
using Base.Test

function test_applying_the_fundamental_solution()
    setprecision(80)do
        info("Testing applying the fundamental solution to the heat equation.")
        #julia> hash(now())
        #     0x9e5efc05faef31aa
        srand(0x9e5efc05faef31aa)
        r5 = region_5(100)
        s5 = [s⁺(r5[i,2],r5[i,3]) for i in 1:size(r5,1)]
        n5 = rand(1:10,size(r5,1))

        r̃₀(s) = M₀(s)*big(e)^(1/(6(abs(s)-0.66)))
        @test all([abs(r₀(n5[i],s5[i])) ≤ abs(r̃₀(s5[i])) for i in 1:size(r5,1)])

    end
end

test_applying_the_fundamental_solution()
