""" s⁺(x,y)

    Return (1-y+im*x)/2 = (1+im*(x+im*y))/2. The superscript + refers to the sign of im*(x+im*y), hence the sign of x in the result.

    This relationship between z=x+im*y and s is taken as standard here since it relates H₀ to Riemann's ξ: H₀(z)=ξ((1+im*z)/2)/8. Thus
    
    (1-y+im*x)/2 = s⁺(x,y)
    (1-y-im*x)/2 = s⁺(x,y)'
    (1+y-im*x)/2 = 1-s⁺(x,y)
    (1+y+im*x)/2 = 1-s⁺(x,y)' (where ' indicates conjugation.)
    """
s⁺(x::Real, y::Real) = (1-y+im*x)/2

s⁺(z::Number) = s⁺(real(z),imag(z))


function M₀(s::Number)
    return 1/8*s*(s-1)/2*π^(-s/2)*√(2*π)*bigexp((s-1)/2*log(s/2)-s/2)
end

function logM₀(s::Number)
    return log(s) + log(s-1) - s/2*log(π) + log(√(2*π)/16) + (s/2-1/2)*log(s/2)-s/2
end

"""
    equation (9), pp 3
    """
function logM₀′(s::Number)
    return 1/(2*s) + 1/(s-1) + 1/2*log(s/(2*π))
end

function logMₑ′(σ::Real, T::Real)
    return logM₀′(σ+im*T)
end

""" α

    alias for logM₀′; equation (9) 
    """
α = logM₀′

function Mₜ(t::Real, s::Number)
    return bigexp(t/4*(α(s))^2)*M₀(s)
end

function Mₜ(t::Real, σ::Real, T::Real)
    return Mₜ(t,σ+im*T)
end

""" B₀(t,x,y)

    Returns Mₜ(t,(1+y-im*x)/2) = Mₜ(t,1-s⁺(x,y)).
    """
function B₀(t::Real, x::Real, y::Real)
    # Recall that (1+y-im*x)/2 = 1-s⁺(x,y).
    return Mₜ(t,1-s⁺(x,y))
end

function B₀(t::Real, z::Number)
    return B₀(t,real(z),imag(z))
end


""" region_5(n)

    Return an nx3 array of n points t,x,y in region 5.
    """
function region_5(n::Integer; xmax::Integer=2000)
    return hcat(rand(n)/2, 200+(xmax-200)*rand(n), rand(n))
end

function in_region_5(t::Real, x::Real, y::Real)
    return (0<t≤1/2) && (x≥200) && (0≤y≤1)
end

function in_region_5(t::Real, z::Number)
    return in_region_5(t, real(z), imag(z))
end

"""
    Equation (14) pp4
    """
function fₜ(t::Real, x::Real, y::Real)
    f1 = f2 = 0.0
    for n in 1:N(t,x)
        b = bᵗₙ(t, n)
        f1 += b/n^(s_star(t,x,y))
        f2 += n^y*b/n^(s_star(t,x,y)'+κ(t,x,y))
    end
    return f1+γₜ(t,x,y)*f2
end

function fₜ(t::Real, z::Number)
    return fₜ(t, real(z), imag(z))
end

"""
    Equation (15) pp 4.
    """
function bᵗₙ(t::Real, n::Int)
    return bigexp(t/4*log(n)^2)
end

"""
    Returns Mₜ((1-y+im*x)/2)/Mₜ((1+y-ix)/2) = Mₜ(s⁺(x,y))/Mₜ(1-s⁺(x,y))
    Equation (16) pp 4

    Named γₜ to avoid naming conflict with Euler's constant
    """
function γₜ(t::Real, x::Real, y::Real)
    # Recall that
    # (1-y+im*x)/2 = s⁺(x,y)
    # (1+y-im*x)/2 = 1-s⁺(x,y)
    s = s⁺(x,y)
    return  Mₜ(t,s)/Mₜ(t,1-s)
end

function γₜ(t::Real, z::Number)
    return γₜ(t, real(z), imag(z))
end

"""
    Equation (17) pp 4.
    """
function s_star(t::Real, x::Real, y::Real)
    # Recall that
    # (1-y+im*x)/2 = s⁺(x,y)
    # (1+y-im*x)/2 = 1-s⁺(x,y)
    s = s⁺(x,y)
    return 1-s + t/2*α(1-s)
end

function s_star(t::Real, z::Number)
    return s_star(t, real(z), imag(z))
end

""" κ(t,x,y)

    Returns t/2*( α((1-y+ix)/2) - α((1+y+ix)/2))
    """
function κ(t::Real, x::Real, y::Real)
    # Recall
    # (1-y+im*x)/2 = s⁺(x,y)
    # (1+y+im*x)/2 = 1-s⁺(x,y)'
    s = s⁺(x,y)
    return t/2*(α(s)-α(1-s'))
end

function κ(t::Real, z::Number)
    return κ(t, real(z), imag(z))
end

function N(t::Real, x::Real)
    return floor(Int, √(x/(4*π) + t/16))
end



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
    return bᵗₙ(t,n)/n^sᵣ*(bigexp((t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66))-1)
end

function bound23(t::Real, x::Real, y::Real)
    bound=0.0
    N₀ = N(t,x)
    γₐ = abs(γₜ(t,x,y))
    κₐ = abs(κ(t,x,y))
    sᵣ = real(s_star(t,x,y))
    for n in 1:N₀
        bound += (1+γₐ*N₀^κₐ*n^y)*bᵗₙ(t,n)/n^sᵣ*ebound_util(n,t,x,y,sᵣ=sᵣ)
    end
    return eA(t,x,y)+eB(t,x,y) ≤ bound
end

function bound23(t::Real, z::Number)
    return bound23(t,real(z),imag(z))
end

function bound24(t::Real, x::Real, y::Real)
    return eC0(t,x,y) ≤ (x/(4*π))^(-(1+y)/4) *
        bigexp(-t/16*log(x/(4*π))^2 + (1.24*(3^y+3^(-y)))/(N(t,x)-0.125) +
               (3*abs(log(x/(4*π))+im*π/2)+10.44)/(x-8.52))
end

function bound24(t::Real, z::Number)
    return bound24(t,real(z),imag(z))
end

