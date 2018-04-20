
using DBNUpperBound
using DBNUpperBound.Asymptotics
using DBNUpperBound.PM15a
using Base.Test

function test_PM15a()

    x = 300.0
    y = 0.4
    z = x+im*y
    s = 1-im*z
    t = 0.4

    @test s = 1+y-im*x
    @test M₀(s) ≈ bigexp(logM₀(s))
    @test logM₀(s)+logM₀′(s)*1e-6 ≈ logM₀(s+1e-6)
    @test logM₀(s)+P.logM₀′(s)*(1e-6)*(1+im) ≈ P.logM₀(s+1e-6*(1+im))
    @test B₀(t, z) == B₀(t, x,y)

    # julia> UInt(Dates.value(Dates.now()))
    # 0x000039e5fab138e8
    srand(0x000039e5fab138e8) #random seed
    r5 = region_5(50)

    all([in_region_5(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    all([bound20(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    all([bound21(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])
    all([bound22(r5[i,1],r5[i,2],r5[i,3]) for i in 1:50])

end
