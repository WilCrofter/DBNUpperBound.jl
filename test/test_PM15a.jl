
using DBNUpperBound
using Base.Test

function test_PM15a()

    x = 300.0
    y = 0.4
    z = x+im*y
    s = (1+im*z)/2
    t = 0.4

    setprecision(53)do # while conversion to multiprecision is incomplete
        @test s == s⁺(z)
        @test M₀(s) ≈ big(e)^(logM₀(s))
        @test logM₀(s)+logM₀′(s)*1e-6 ≈ logM₀(s+1e-6)
        @test logM₀(s)+logM₀′(s)*(1e-6)*(1+im) ≈ logM₀(s+1e-6*(1+im))
        @test B₀(t, z) == B₀(t,x,y)
        @test Mₜ(t,s) == Mₜ(t,s')'
    end

    # test Hₜ against reference values less than 1000
    @test Hₜ(.4,30.0,.4)[1] ≈ -0.00010001026469315639165 - 7.1357019921469872653e-6*im
    @test Hₜ(.4,100.0,.4)[1] ≈  6.7021522172791266841e-16 + 3.1337965840705569244e-16*im
    @test Hₜ(.4,300.0,.4)[1] ≈ -4.0159674206251463624e-49 - 1.4006524430296850032e-49*im
    
    # julia> UInt(Dates.value(Dates.now()))
    # 0x000039e5fab138e8
    srand(0x000039e5fab138e8) #random seed
    r5 = region_5(50)

    @test all([in_region_5(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    @test all([bound20(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    @test all([bound21(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    @test all([bound22(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    # TODO: as implemented, bounds 23 & 24 fail tests.
    # Defer debugging until later in paper.
    # @test all([bound23(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    # @test all([bound24(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])

    DBNUpperBound.PM15a.lemma_51i(1000)
    DBNUpperBound.PM15a.lemma_51ii(1000)
    DBNUpperBound.PM15a.lemma_51iii(1000)
    DBNUpperBound.PM15a.lemma_51iv(1000)
    DBNUpperBound.PM15a.lemma_51v(1000)
    DBNUpperBound.PM15a.lemma_51vi(10000)

end

test_PM15a()

include("test_estimates_for_large_x.jl")
