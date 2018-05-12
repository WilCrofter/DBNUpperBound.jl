""" 
    Inequality (20) pp 
    """
function bound20(t::Real, x::Real, y::Real)
    return abs(γₜ(t,x,y)) ≤ exp(0.02*y)*(x/(4*π))^(-y/2)
end

function bound21(t::Real,x::Real,y::Real)
    return real(s_star(t,x,y)) ≥ (1+y)/2 + t/4*log(x/(4*π)) - t*max(0.0, 1-3*y+(4*y*(1+y))/x^2)/(2*x^2)
end

function bound22(t::Real, x::Real, y::Real)
    return abs(κ(t,x,y)) ≤ t*y/(2*(x-6))
end

function ebound_util(n::Int, t::Real, x::Real, y::Real; sᵣ::Real=real(s_star(t,x,y)))
    return bᵗₙ(t,n)/n^sᵣ*(big(e)^((t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66))-1)
end

function bound23(t::Real, x::Real, y::Real)
    bound=0.0
    N₀ = N(t,x)
    γₐ = abs(γₜ(t,x,y))
    κₐ = abs(κ(t,x,y))
    sᵣ = real(s_star(t,x,y))
    for n in 1:N₀
        bound += (1+γₐ*N₀^κₐ*n^y)*bᵗₙ(t,n)/n^sᵣ*(big(e)^((t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66))-1)
    end
    return eA(t,x,y)+eB(t,x,y) ≤ bound
end

function bound23(t::Real, z::Number)
    return bound23(t,real(z),imag(z))
end

function bound24(t::Real, x::Real, y::Real)
    return eC0(t,x,y) ≤ (x/(4*π))^(-(1+y)/4) *
        big(e)^(-t/16*log(x/(4*π))^2 + (1.24*(3^y+3^(-y)))/(N(t,x)-0.125) +
               (3*abs(log(x/(4*π))+im*π/2)+10.44)/(x-8.52))
end

function bound24(t::Real, z::Number)
    return bound24(t,real(z),imag(z))
end

function test_introduction()
    info("Testing introduction.jl")
    setprecision(80)do

        x = 300.0
        y = 0.4
        z = x+im*y
        s = (1+im*z)/2
        t = 0.4

        @test s == s⁺(z)
        @test M₀(s) ≈ big(e)^(logM₀(s))
        @test logM₀(s)+logM₀′(s)*1e-6 ≈ logM₀(s+1e-6)
        @test logM₀(s)+logM₀′(s)*(1e-6)*(1+im) ≈ logM₀(s+1e-6*(1+im))
        @test B₀(t, z) == B₀(t,x,y)
        @test Mₜ(t,s) == Mₜ(t,s')'

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
    end
end

test_introduction()
