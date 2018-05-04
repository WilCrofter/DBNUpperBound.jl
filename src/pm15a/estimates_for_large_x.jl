export ϵₜₙ, rₜₙ, ϵ̃, eA, eB, eC, eC0

#= Recall
    (1-y+im*x)/2 = s⁺(x,y)
    (1-y-im*x)/2 = s⁺(x,y)'
    (1+y-im*x)/2 = 1-s⁺(x,y)
    (1+y+im*x)/2 = 1-s⁺(x,y)'
=#

""" 
    Definition (43) pp 14
    """
function ϵₜₙ(t::Real, n::Int, σ::Real, T::Real)
    return bigexp((t^2/8*abs(α(σ+im*T)-log(n))^2+t/4+1/6)/(T-3.33))-1
end

function ϵₜₙ(t::Real, n::Int, u::Number)
    return ϵₜₙ(t, n, real(u), imag(u))
end

""" rₜₙ(t::Real, n::Int, σ::Real, T::Real)

    Estimate of rₜₙ, Proposition 6.1 pp 13
    Returns the basic estimate and error term ϵₜₙ(σ+iT).
    """
function rₜₙ(t::Real, n::Int, σ::Real, T::Real)
    0 < t ≤ 1/2 || error("t must be positive and at most 1/2")
    10 < T || error("T must be at least 10.")
    s=σ+im*T
    r = Mₜ(t,σ,T)*bᵗₙ(t, n)*bigexp(-log(n)*(s+t/2*α(s)))
    return r, ϵₜₙ(t,n,σ,T)
end                   

function rₜₙ(t::Real, n::Int, s::Number)
    return rₜₙ(t, n, real(s), imag(s))
end

"""
    Definition (51) pp 15
    """
function C₀(p::Real)
    return p ≈ 0.5 || p ≈ -0.5 ? (1-im)/4 :(exp(π*im*(p^2/2 + 3/8)) - im*√(2)*cos(π*p/2))/2*cos(π*p)
end


"""
    Definition (57) pp 16
    """
function ϵ̃(t::Real, σ::Real, T::Real)
    T′ = T - π*t/8
    a = √(T′/(2*π))
    return ((0.397*9^σ)/(a-0.125) + 5/(3*(T′-3.33)))*bigexp(3.49/(T′-3.33))
end

function ϵ̃(t::Real, u::Number)
    return ϵ̃(t,real(u),imag(u))
end


"""
    Definition (69) pp. 23.
    """
function eA(t::Real, x::Real, y::Real)
    #= Recall
    (1-y+im*x)/2 = s⁺(x,y)
    =#
    s = s⁺(x,y)
    ans = 0.0
    skᵣ = real(s)+real(κ(t,x,y))
    for n in 1:N(t,x)
        ans += n^y*(bᵗₙ(t,n)/n^skᵣ)*ϵₜₙ(t,n,s)
    end
    return abs(γₜ(t,x,y))*ans
end

function eA(t::Real, z::Number)
    return eA(t, real(z), imag(z))
end

"""
    Definition (70) pp. 23.
    """
function eB(t::Real, x::Real, y::Real)
    #= Recall
    (1+y+im*x)/2 = 1-s⁺(x,y)'   
    =#
    s= s⁺(x,y)
    ans = 0.0
    for n in 1:N(t,x)
        ans += bᵗₙ(t,n)/n^(real(s))*ϵₜₙ(t,n,1-s')
    end
    return ans
end

function eB(t::Real, z::Number)
    return eB(t,real(z),imag(z))
end

"""
    Definition (71) pp. 23.
    """
function eC(t::Real, x::Real, y::Real)
    #= Recall
    (1-y+im*x)/2 = s⁺(x,y)
    (1+y+im*x)/2 = 1-s⁺(x,y)'
    =#
    s = s⁺(x,y)
    T = imag(s)
    T′ = T - π*t/8
    return exp(t*π^2/64)*abs(M₀(im*T′))/abs(Mₜ(t, 1-s'))*(ϵ̃(t,s) + ϵ̃(t,1-s'))
end

function eC(t::Real, z::Number)
    return eC(t,real(z),imag(z))
end


"""
    Definition (72) pp. 23.
    """
function eC0(t::Real, x::Real, y::Real)
    #= Recall
    (1-y+im*x)/2 = s⁺(x,y)
    (1+y+im*x)/2 = 1-s⁺(x,y)'
    =#
    s = s⁺(x,y)
    T = imag(s)
    T′ = T - π*t/8
    return exp(t*π^2/64)*abs(M₀(im*T′))/abs(Mₜ(t,1-s'))*(1 + ϵ̃(t,s) + ϵ̃(t,1-s'))
end

function eC0(t::Real, z::Number)
    return eC0(t,real(z),imag(z))
end

