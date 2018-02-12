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
    function fKKL(u) u^(1/4)*(2*ψKKL(u)+1) end
    @test all([fKKL(u) ≈ fKKL(1/u) for u in x])

    # Verify that ϕKKL(x) ≈ ϕKKL(1/x) unless their difference is less than 1e-14
    # which seems the best we can do with Float64.
    @test all([isapprox(ϕKKL(u),ϕKKL(1/u),atol=1e-14) for u in x])

    # Verify that ϕKKL(big(x)) ≈ ϕKKL(1/big(x)) (multiprecision)
    @test all([isapprox(ϕKKL(big(u)),ϕKKL(1/big(u)),atol=1e-16) for u in x])

    # Verify that two implementations of ϕKKL are approximately equal unless
    # their difference is less than 1e-14, which seems the best we can
    # do with Float64. As above, multiprecision does slightly better but not much.
    @test all([isapprox(ϕKKL(u),DBNUpperBound.ϕKKL4test(u),atol=1e-14) for u in x])

    # KKL2009 claims that ΞKKL(0.0,t) = ΞRL(t) = ξ(1/2 + t*im), but as those functions
    # are implemented here, following KKL2009 for implementation of ΞKKL(0.0,t)
    # and using a standard definition of Riemann's xi for implementation
    # of ξ, the actual relation seems to be ΞKKL(0.0,t) = ξ(1/2 + t*im)/4 = ΞRL(t)/4.
    # where  ΞRL(t) is the Riemann-Landau Xi.
    @test all([ΞKKL(0.0,t)[1] ≈ ξ(1/2+im*t)/4 for t in x])
    # Equivalently
    @test all([ΞKKL(0.0,t)[1] ≈ ΞRL(t)/4 for t in x])

    # Assuming ΞKKL(0.0,t) = ξ(1/2 + t*im)/4 where Ξ(0.0,t) is as defined in KKL2009
    # and noting that, as Ht(0.0,t) is defined at PM15, Ht(0.0,t)= H0(t) =ξ((1+t*im)/2)/8,
    # and that ξ((1+t*im)/2)/8 = 1/2*ξ(1/2+t/2*im))/4 = 1/2*Ξ(0.0,t/2), we should
    # have Ht(0.0,t)=H0(t)=1/2*ΞKKL(0.0,t/2). We know from test_Ht.jl that
    # the first equality holds within numerical tolerance. Checking the second equality:
    @test all([H0(t)≈1/2*ΞKKL(0.0,t/2)[1] for t in x])

    # Using Ht(0.0,t)=1/2*ΞKKL(0.0,t/2), change of variable, and Fourier inversion
    # we can deduce that Φpm(t)≈1/2*ΦKKL(2*t), where Φpm is Φ as defined in PM15,
    # and ΦKKL as defined in KKL2009
    all([Φpm(t)≈1/2*ΦKKL(2*t) for t in x])
end

test_KKL()
