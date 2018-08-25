export s⁺, M₀, logM₀, logM₀′, α, α′, Mₜ, B₀, bᵗₙ, logbᵗₙ
export region_5, in_region_5
export γₜ, κ, fₜ
export N, H̃, H̃₂

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

"""

    M₀(s) := 1/8 × s(s-1)/2 × π^(-s/2)√(2π) × exp((s/2-1/2)Log(s/2)-s/2) Equation (6) pp 18.
    """
function M₀(s::Number)
    return big(e)^(logM₀(s))
end

function logM₀(s::Number)
    if s ≈ 0 || s ≈ 1
        return -Inf
    else
        return log(s) + log(s-1) - s/2*log(π) + log(√(2*π)/16) + (s/2-1/2)*log(s/2)-s/2
    end
end

"""
    equation (9), pp 3
    """
function logM₀′(s::Number)
    return 1/(2*s) + 1/(s-1) + 1/2*log(s/(2*π))
end

function logM₀′(σ::Real, T::Real)
    return logM₀′(σ+im*T)
end

""" α

    alias for logM₀′; equation (9) 
    """
α = logM₀′

""" α′

    Derivative of α
"""
α′(s::Number) = -1/(2*s^2)-1/(s-1)^2 + 1/(2*s)

function Mₜ(t::Real, s::Number)
    return big(e)^(t/4*(α(s))^2)*M₀(s)
end

function Mₜ(t::Real, σ::Real, T::Real)
    return Mₜ(t,σ+im*T)
end

""" B₀(t,x,y)

    Returns Mₜ(t,(1+y-im*x)/2) = Mₜ(t,1-s⁺(x,y)) as in Def (11), pp. 3 
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
    return hcat(rand(n)/2, 200 .+ (xmax-200)*rand(n), rand(n))
end

function in_region_5(t::Real, x::Real, y::Real)
    return (0≤t≤1/2) && (x≥200) && (0≤y≤1)
end

function in_region_5(t::Real, z::Number)
    return in_region_5(t, real(z), imag(z))
end

""" fₜ(t::Real, x::Real, y::Real)

    Returns fₜ, computed as follows. We have fₜ = A/B₀ + B/B₀ by (14) pp. 4, the definitions of A and B in Cor 6.4 pp. 22, the definition of B₀ in (11) pp. 3, and the definition of γ in (16) pp. 4. 
    """
function fₜ(t::Real, x::Real, y::Real)
    ans = H̃(t,x,y)
    B0 = B₀(t,x,y)
    return ans[1]/B0, ans[2]/abs(B0)
end

function fₜ(t::Real, z::Number)
    return fₜ(t, real(z), imag(z))
end

"""
    Equation (15) pp 4.
    """
function bᵗₙ(t::Real, n::Int)
    return big(e)^(t/4*log(n)^2)
end

function logbᵗₙ(t::Real, n::Int)
    return t/4*log(n)^2
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

"""  H̃(t::Real, x::Real, y::Real)
    
    Returns the "A+B" approximation to Hₜ and the error bound, EA+EB+EC0. 
    """
function H̃(t::Real, x::Real, y::Real)
    Ã,EA = A(t,x,y)
    B̃,EB = B(t,x,y)
    EC0  = C(t,x,y)[3]
    return Ã+B̃, EA+EB+EC0
end

""" H̃₂(t::Real, x::Real, y::Real)
    
    Returns the "A+B-C" approximation to Hₜ and the error bound, EA+EB+EC
    """
function H̃₂(t::Real, x::Real, y::Real)
    Ã,EA = A(t,x,y)
    B̃,EB = B(t,x,y)
    C̃,EC,EC0  = C(t,x,y)
    return Ã+B̃-C̃, EA+EB+EC
end
 

