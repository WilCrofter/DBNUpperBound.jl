
using DBNUpperBound
using QuadGK
using Base.MathConstants
using SpecialFunctions
using Test
using Random

function test_estimating_sums()
    @info("Testing estimating sums")
    setprecision(80) do
        N₀ = 10^4 # small to test times low
        y = t = .2
        # N = √(x/(4*π) + t/16)
        X = 4*π*(N₀^2-t/16)
        w = ζ = 1.0
        fn₀ = FN(X + im*y, ζ, w, N₀)
        N₁ = 2*N₀ # to keep test times low
        H = round(Int,(N₁-N₀)/3)    # ditto
        T = 10
        centers = NᵢRange(N₁,N₀,H)
        # compute approximate partial sums 1 at a time
        nsums = length(centers)
        psums = Vector{Complex{BigFloat}}(undef, nsums)
        errs = Vector{BigFloat}(undef, nsums)
        fullsums = Vector{Complex{BigFloat}}(undef, nsums)
        dN = maximum(DBNUpperBound.PM15a.hRange(H)) # hRange is not exported
        for i in 1:nsums
            psums[i],errs[i] = FN_tail(X+im*y, ζ, w, N₁, N₀, H, T; centers=centers[i:i])
            fullsums[i] = FN(X+im*y,ζ,w,centers[i].+ dN)
        end
        psums = fn₀ .+ cumsum(psums)
        errs = cumsum(errs)
        # test that exact sums are bounded in absolute value by approximations + errors
        for i in 1:nsums
            @test abs(fullsums[i]) ≤ abs(psums[i])+errs[i]
            @test abs(fullsums[i]-psums[i]) ≤ errs[i]
        end
    end
end


test_estimating_sums()
