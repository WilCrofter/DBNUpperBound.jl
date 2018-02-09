using DBNUpperBound
using QuadGK
using SpecialFunctions
using Base.Test

function test_KKL()
    # julia> rand(UInt32) # when test was written
    # 0xdf72f388
    srand(0xdf72f388)
    x = rand(100) # 100 random numbers btw 0 and 1

    # verify functional equation for ψ
    function fKKL(u) u^(1/4)*(2*ψ(u)+1) end
    @test all([fKKL(u) ≈ fKKL(1/u) for u in x])

    # Verify that ϕ(x) ≈ ϕ(1/x) unless their difference is less than 1e-14
    # which seems the best we can do with Float64.
    @test all([isapprox(ϕ(u),ϕ(1/u),atol=1e-14) for u in x])

    # Verify that ϕ(big(x)) ≈ ϕ(1/big(x)) (multiprecision)
    @test all([isapprox(ϕ(big(u)),ϕ(1/big(u)),atol=1e-16) for u in x])

    # Verify that two implementations of ϕ are approximately equal unless
    # their difference is less than 1e-14, which seems the best we can
    # do with Float64. As above, multiprecision does slightly better but not much.
    @test all([isapprox(ϕ(u),DBNUpperBound.ϕ4test(u),atol=1e-14) for u in x])

    # KKL2009 claims that Ξ(0.0,t) = ξ(1/2 + t*im), but as those functions
    # are implemented here, following KKL2009 for implementation of Ξ(0.0,t)
    # and the using a standard definition of Riemann's xi, for implementation
    # of ξ, the actual relation seems to be Ξ(0.0,t) = ξ(1/2 + t*im)/4
    @test all([Ξ(0.0,t)[1] ≈ ξ(1/2+im*t)/4 for t in x])
end

test_KKL()
