
using DBNUpperBound
using DBNUpperBound.Asymptotics
using QuadGK
using SpecialFunctions
using Base.Test


function test_asymptotics()

    t = 0.4
    y = 0.4
    x = [10.0^n for n in 5:2:15] # For n much bigger than 15, Inf's will appear.
    z = x+im*y
    s = (1+im*z)/2

    function ABC(s)
        M=N=floor(Int,sqrt(imag(s)/(2*π)))
        a=A(t,N,s)
        b=B(t,M,s)
        c=C(t,N,M,s)
        return a,b,c
    end

    abc=[ABC(u) for u in s]

    # Test that B(t,M,s)/B0 → 1 as x→∞ where B0=B(t,1,s)
    # or, actually, that the ratios rounded to the nearest integer are 1 for
    # reasonably large x. Convergence apparently isn't monotone or particularly fast.

    for i in eachindex(abc)
        if x[i] > x[2]
            b = abc[i][2]
            @test round(abs(b/B(t,1,s[i])))==1
        end
    end

end

test_asymptotics()
