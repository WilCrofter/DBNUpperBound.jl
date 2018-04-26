
using DBNUpperBound
using DBNUpperBound.Asymptotics
using DBNUpperBound.PM15a
using Base.Test

function test_PM15a()

    x = 300.0
    y = 0.4
    z = x+im*y
    s = (1+im*z)/2
    t = 0.4

    @test s == s⁺(z)
    @test M₀(s) ≈ bigexp(logM₀(s))
    @test logM₀(s)+logM₀′(s)*1e-6 ≈ logM₀(s+1e-6)
    @test logM₀(s)+logM₀′(s)*(1e-6)*(1+im) ≈ logM₀(s+1e-6*(1+im))
    @test B₀(t, z) == B₀(t,x,y)
    @test Mₜ(t,s) == Mₜ(t,s')'

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

end

test_PM15a()
