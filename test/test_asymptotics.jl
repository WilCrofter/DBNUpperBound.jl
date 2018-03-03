
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

    #= Branch or mod 2π problem?
    for u in s
        α = u/(2*π*im*Ns(s))-Ns(s)
        # \Psi(\alpha) = 2\pi i - e^{-2\pi i \alpha} \Psi(\alpha+1)
        @test Ψpm(α) ≈ 2*π*im - exp(-2*π*im*α)*Ψpm(α+1)
        # \Psi(\alpha+1)-\Psi(\alpha) = \frac{2\pi}{\sqrt{i}} \exp(i \pi \alpha^2)
        @test Ψpm(α+1)-Ψpm(α) ≈ (2*π/sqrt(im))*exp(im*π*α^2)
    end
    =#

    # Test B0, B0', B0eff against data published at Polymath15 test problem:
    # http://michaelnielsen.org/polymath1/index.php?title=Polymath15_test_problem
    data=((1e3,[(3.4405+3.5443im)*big(1)e-167
                (3.4204+3.5383im)*big(1)e-167
                (3.4426+3.5411im)*big(1)e-167
                (2.3040+2.3606im)*big(1)e-167]),
          (1e4,[(-1.1843-7.7882im)*big(1)e-1700
	        (-1.1180-7.7888im)*big(1)e-1700
	        (-1.1185-7.7879im)*big(1)e-1700
	        (-1.1155-7.5753im)*big(1)e-1700]),
          (1e5,[(-7.6133+2.5065im)*big(1)e-17047
	        (-7.6134+2.5060im)*big(1)e-17047
	        (-7.6134+2.5059im)*big(1)e-17047
	        (-7.5483+2.4848im)*big(1)e-17047]),
          (1e6,[(-3.1615-7.7093im)*big(1)e-170537
	        (-3.1676-7.7063im)*big(1)e-170537
	        (-3.1646-7.7079im)*big(1)e-170537
	        (-3.1590-7.6898im)*big(1)e-170537]),
          (1e7,[(2.1676-9.6330im)*big(1)e-1705458
	        (2.1711-9.6236im)*big(1)e-1705458
	        (2.1571-9.6329im)*big(1)e-1705458
	        (2.2566-9.6000im)*big(1)e-1705458]))
    for datum in data
        x = datum[1]
        z = x+im*.4
        s = (1+im*z)/2
        @test isapprox(B(.4,1,s),datum[2][1],rtol=4)
        @test isapprox(Bprime(.4,s,1),datum[2][2],rtol=4)
        @test isapprox(Beff(.4,s,1),datum[2][3],rtol=4)
    end
end

test_asymptotics()

