function M₀(s::Number)
    return 1/8*s*(s-1)/2*π^(-s/2)*√(2*π)*bigexp((s-1)/2*log(s/2)-s/2)
end

function logM₀(s::Number)
    return log(s) + log(s-1) - s/2*log(π) + log(√(2*π)/16) + (s/2-1/2)*log(s/2)-s/2
end

function logM₀′(s::Number)
    return 1/(2*s) + 1/(s-1) + 1/2*log(s/(2*π))
end

""" α

    alias for logM₀′
    """
α = logM₀′

function Mₜ(t::Real, s::Number)
    return bigexp(t/4*(logM₀′(s))^2)*M₀(s)
end

function B₀(t::Real, z::Number)
    return Mₜ(t,(1-im*z)/2)
end

function B₀(t::Real, x::Real, y::Real)
    return Mₜ(t,(1+y-im*x)/2)
end

""" region_5(n)

    Return an nx3 array of n points t,x,y in region 5.
    """
function region_5(n::Integer; Nmax::Integer=2000)
    return hcat(rand(n)/2, 200+(Nmax-200)*rand(n), rand(n))
end

function in_region_5(t::Real, x::Real, y::Real)
    return (0<t≤1/2) && (x≥200) && (0≤y≤1)
end

function in_region_5(t::Real, z::Number)
    return in_region_5(t, real(z), imag(z))
end

function bᵗₙ(t::Real, n::Int)
    return bigexp(t/4*log(n)^2)
end

function γ(t::Real, x::Real, y::Real)
    # Definition in the Apr 19 version of the paper
    # is Mₜ(t,(1-y+im*x)/2)/Mₜ(t,(1+y-im*x)/2)
    # which seems to be a typo. See proof of Prop 6.6.
    return  Mₜ(t,(1-y+im*x)/2)/Mₜ(t,(1+y+im*x)/2)
end

function γ(t::Real, z::Number)
    return γ(t, real(z), imag(z))
end

function s_star(t::Real, x::Real, y::Real)
    return (1+y-im*x)/2 + t/2*α((1+y-im*x)/2)
end

function s_star(t::Real, z::Number)
    return s_star(t, real(z), imag(z))
end

function κ(t::Real, x::Real, y::Real)
    return t/2*(α((1-y+im*x)/2)-α((1+y+ix)/2))
end

function κ(t::Real, z::Number)
    return κ(t, real(z), imag(z))
end

function N(t::Real, x::Real)
    return floor(Int, √(x/(4*π) + t/16))
end

function fₜ(t::Real, x::Real, y::Real)
    f1 = f2 = 0.0
    for n in 1:N(t,x)
        f1 += bᵗₙ(t, n)
        f2 += n^y*bᵗₙ(t, n)/n^(s⋆(t,x,y)'+κ(t,x,y))
    end
    return f1+γ(t,x,y)*f2
end

function fₜ(t::Real, z::Number)
    return fₜ(t, real(z), imag(z))
end

function bound20(t::Real, x::Real, y::Real)
    return abs(γ(t,x,y)) ≤ exp(0.02*y)*(x/(4*π))^(-y/2)
end

function bound21(t::Real,x::Real,y::Real)
    return real(s_star(t,x,y)) ≥ (1+y)/2 + t/4*log(x/(4*π)) - t*max(0.0, 1-3*y+(4*y*(1+y))/x^2)
end

function bound22(t::Real, x::Real, y::Real)
    return abs(κ(t,x,y)) ≤ t*y/(2*(x-6))
end

function ϵ̃(t::Real, σ::Real, T::Real)
    T′ = T - π*t/8
    a = √(T′/(2*π))
    return ((0.397*9^σ)/(a-0.125) + 5/(3*(T′-3.33)))*bigexp(3.49/(T′-3.33))
end

function ϵ̃(t::Real, u::Number)
    return ϵ̃(t,real(u),imag(u))
end

function ϵₜₙ(t::Real, σ::Real, T::Real)
    return bigexp((t^2/8*abs(α(σ+im*T)-log(n))^2+t/4+1/6)/(T-3.33))-1
end

function ϵₜₙ(t::Real, u::Number)
    return ϵₜₙ(t, real(u), imag(u))
end
