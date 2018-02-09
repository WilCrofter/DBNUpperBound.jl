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

    # KKL2009 claims that Ξ(0.0,t) = Ξ(t) = ξ(1/2 + t*im), but as those functions
    # are implemented here, following KKL2009 for implementation of Ξ(0.0,t)
    # and using a standard definition of Riemann's xi for implementation
    # of ξ, the actual relation seems to be Ξ(0.0,t) = ξ(1/2 + t*im)/4 = Ξ(t)/4.
    # where  Ξ(t) is the Riemann-Landau Xi.
    @test all([Ξ(0.0,t)[1] ≈ ξ(1/2+im*t)/4 for t in x])
    # Equivalently
    @test all([Ξ(0.0,t)[1] ≈ Ξ(t)/4 for t in x])

    # Assuming Ξ(0.0,t) = ξ(1/2 + t*im)/4 where Ξ(0.0,t) is as defined in KKL2009
    # and noting that, as Ht(0.0,t) is defined at PM15, Ht(0.0,t)= H0(t) =ξ((1+t*im)/2)/8,
    # and that ξ((1+t*im)/2)/8 = 1/2*ξ(1/2+t/2*im))/4 = 1/2*Ξ(0.0,t/2), we should
    # have Ht(0.0,t)=H0(t)=1/2*Ξ(0.0,t/2). We know from test_Ht.jl that
    # the first equality holds within numerical tolerance. Checking the second equality:
    @test all([H0(t)≈1/2*Ξ(0.0,t/2)[1] for t in x])
end

test_KKL()
