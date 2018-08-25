
export ϵₜₙ, ϵ̃ₜₙ, ϵ̃, rₜₙ, RtN, A, B, C, EA, EB, EC, EC₀
export eA, eB, eC, eC₀
export Aterm, Bterm

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
    return big(e)^((t^2/8*abs(α(σ+im*T)-log(n))^2+t/4+1/6)/(T-3.33))-1
end

function ϵₜₙ(t::Real, n::Int, u::Number)
    return ϵₜₙ(t, n, real(u), imag(u))
end

"""
    Upper bound for ϵₜₙ as given in the last line of the proof of Prop 6.6 iv.
    ϵₜₙ((1±y+ix)/2) ≤ exp((t²/32*log²(x/(4πn²))+0.313)/(T-3.33))-1
    equivalently, ϵₜₙ((1±y+ix)/2) ≤ exp((t²/16*log²(x/(4πn²))+0.626)/(x-6.66))-1
    since T=x/2 by convention.

    Note that x-6.66 (as in 6.6 v) rather than x/2-6.66 (as in 6.6 iv) appears in the denominator here. The x/2 is assumed to be a typo meant as T/2 (=x.)
    """
function ϵ̃ₜₙ(t::Real, n::Int, σ::Real, T::Real)
    x=2*T
    return big(e)^((t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66)) - 1
end

function ϵ̃ₜₙ(t::Real, n::Int, s::Number)
    return  ϵ̃ₜₙ(t,n,real(s),imag(s))
end

""" rₜₙ(t::Real, n::Int, σ::Real, T::Real)

    Estimate of rₜₙ, Proposition 6.1 pp 13
    Returns the basic estimate and error term ϵₜₙ(σ+iT).
    """
function rₜₙ(t::Real, n::Int, σ::Real, T::Real)
    0 < t ≤ 1/2 || error("t must be positive and at most 1/2")
    10 < T || error("T must be at least 10.")
    s=σ+im*T
    r = Mₜ(t,s)*bᵗₙ(t, n)*big(e)^(-log(n)*(s+t/2*α(s)))
    return r, 1+ϵₜₙ(t,n,s)
end                   

function rₜₙ(t::Real, n::Int, s::Number)
    return rₜₙ(t, n, real(s), imag(s))
end

function est_utility(t::Real, σ::Real, T::Real)
    T′ = T+π*t/8  # = x/2 +πt/8 assuming σ+iT = (1+i(x+iy))/2 (66) pp 22
    a = √(T′/(2*π)) # (47) pp 15
    # Compare N(t,x) = floor(Int, √(x/(4*π) + t/16)) introduction.jl and (48) pp 15
    # with N₀ = floor(Int,a)
    # The two expressions are idential provided T = x/2 which is the case by convention.
    N₀ = N(t,2*T)
    p = 1-2*(a-N₀) # (49) pp 15
    U = big(e)^(-im*(T′/2*log(T′/(2*π))-T′/2-π/8)) # (50) pp 15
    return T′, a, N₀, p, U
end

function est_utility(t::Real, s::Number)
    return est_utility(t,real(s),imag(s))
end

"""
    Estimate of RₜN Proposition 6.3 pp 16
    """
function RtN(t::Real, σ::Real, T::Real)
    T ≥ 100.0 || error("T = ℑ(s) must be at least 100.0")
    T′, a, N₀, p, U = est_utility(t,σ,T)
    return (-1)^(N₀-1)*U*big(e)^(im*π*σ/4+t*π^2/64)*M₀(im*T′)*(C₀(p)+ϵ̃(σ+im*T))
end

"""
    Definition (51) pp 15
    """
function C₀(p::Real)
    return p ≈ 0.5 || p ≈ -0.5 ? (1-im)/4 : (exp(π*im*(p^2/2 + 3/8)) - im*√(2)*cos(π*p/2))/2*cos(π*p)
end

""" Aterm(t::Real, n::Int, x::Real, y::Real; s=s⁺(x,y), astar=s+t/2*α(s))

    A utility which returns a term bᵗₙ/n^(s+t⋅α(s)/2), useful for computing A and related quantities. Recall that (1-y+im*x)/2 = s⁺(x,y).
    """
function Aterm(t::Real, n::Int, x::Real, y::Real; s::Number=s⁺(x,y), astar::Number=s+t/2*α(s))
    return bᵗₙ(t,n)*big(e)^(-log(n)*astar)
end


""" A(t::Real, x::Real, y::Real)

    Returns A and error estimate EA for A
    where A(x+iy) = Mₜ((1-y+ix)/2)*∑bᵗₙ/n^((1-y+ix)/2 +t/2*α((1-y+ix)/2))
          EA(x+iy) = |Mₜ((1-y+ix)/2)|*∑bᵗₙ/n^((1-y+ix)/2 +t/2*α((1-y+ix)/2))|*ϵₜₙ(t,x,y)
    Corollary 6.4, pp 22. 
    """
function A(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    astar = (s+t/2*α(s))
    ans = err = 0.0
    for n in 1:N(t,x)
        term = Aterm(t,n,x,y,s=s,astar=astar)
        ans += term
        err += abs(term)*ϵₜₙ(t,n,s)
    end
    return Mₜ(t,s)*ans, abs(Mₜ(t,s))*err
end

function A(t::Real, z::Number)
    return A(t, real(z), imag(z))
end

""" Bterm(t::Real, n::Int, x::Real, y::Real; s=s⁺(x,y), bstar=1-s+t/2*α(1-s))

    A utility which returns a term bᵗₙ/n^(1-s+t⋅α(1-s)/2), useful for computing B and related quantities. Recall that (1+y-im*x)/2 = 1-s⁺(x,y).
    """
function Bterm(t::Real, n::Int, x::Real, y::Real; s=s⁺(x,y), bstar=1-s+t/2*α(1-s))
    return bᵗₙ(t,n)*big(e)^(-log(n)*bstar)
end

""" Returns B and error estimate EB for B

    where B(x+iy) = Mₜ((1+y-ix)/2)*∑bᵗₙ/n^((1+y-ix)/2 +t/2*α((1+y-ix)/2))
    and  EB(x+iy) = |Mₜ((1+y-ix)/2)|*∑|bᵗₙ/n^((1+y-ix)/2 +t/2*α((1+y-ix)/2))|*ϵₜₙ((1+y+ix)/2)
    Corollary 6.4, pp 22. 
    """
function B(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    ans = err = 0.0
    # Note that (1+y-ix)/2 = 1-s⁺(x,y) and (1+y+ix)/2 = 1-s⁺(x,y)' where ' indicates conjugation.
    bstar = 1-s + t/2*α(1-s)
    for n in 1:N(t,x)
        term = Bterm(t,n,x,y,s=s,bstar=bstar)
        ans += term
        err += abs(term)*ϵₜₙ(t,n,1-s')
    end
    return Mₜ(t,1-s)*ans, abs(Mₜ(t,1-s))*err
end

""" 
    Returns C and error estimates EC, and EC₀
    where C = 2*(-1)^N*exp(-i*π*y/8 +t*π^2/64)*real(M₀(im*T′)*C₀(p)*U*exp(π*i/8))
         EC = exp(t*π^2/64)*abs(M₀(im*T′))*(ϵ̃(t,s)+ϵ̃(t,1-s'))
    and EC₀ = exp(t*π^2/64)*abs(M₀(im*T′))*(1+ϵ̃(t,s)+ϵ̃(t,1-s'))
    Corollary 6.4 pp. 22
    """
function C(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′ = T+π*t/8  # = x/2 +πt/8 assuming σ+iT = (1+i(x+iy))/2 (66) pp 22
    a = √(T′/(2*π)) # (47) pp 15
    # Compare N(t,x) = floor(Int, √(x/(4*π) + t/16)) introduction.jl and (48) pp 15
    # with N₀ = floor(Int,a)
    # The two expressions are idential provided T = x/2 which is the case by convention.
    N₀ = N(t,2*T)
    p = 1-2*(a-N₀) # (49) pp 15
    U = big(e)^(-im*(T′/2*log(T′/(2*π))-T′/2-π/8)) # (50) pp 15
    ans = 2*(-1)^N₀*big(e)^(-im*π*y/8 +t*π^2/64)*real(M₀(im*T′)*C₀(p)*U*exp(π*im/8))
    tmp = big(e)^(t*π^2/64)*abs(M₀(im*T′))
    ẽ=ϵ̃(t,s)+ϵ̃(t,1-s')
    EC = tmp*ẽ
    EC₀= tmp*(1+ẽ)
    return ans, EC, EC₀
end

""" EA(t::Real,x::Real,y::Real)

    EA(x+iy) :=  |Mₜ((1-y+ix)/2)|*∑bᵗₙ/n^((1-y)/2 +t/2*ℜ(α((1-y+ix)/2)))*ϵₜₙ((1-y+ix)/2)
    Corollary 6.4 pp 22
    """
function EA(t::Real,x::Real,y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    # rastar is the real part of astar = s+t/2*α(s) = (1-y)/2 + t/2*α((1-y+ix)/2)
    # which appears in the definition of A
    rastar = real(s+t/2*α(s))
    σ, T = real(s), imag(s)
    T′, a, N₀, p, U = est_utility(t,σ,T)
    # Recall that (1-y+im*x)/2 = s⁺(x,y)
    ans = big(0.0)
    for n in 1:N₀
        ans += bᵗₙ(t,n)*big(e)^(-log(n)*rastar)*ϵₜₙ(t,n,s)
    end
    return abs(Mₜ(t,s))*ans
end

""" eA(t::Real,x::Real,y::Real)

    Returns EA/|B₀|
    """
eA(t::Real,x::Real,y::Real) = EA(t,x,y)/abs(B₀(t,x,y))

""" EB(t::Real,x::Real,y::Real)

    EB(x+iy) :=  |Mₜ((1+y+ix)/2)|*∑bᵗₙ/n^((1+y)/2 +t/2*ℜ(α((1+y+ix)/2)))*ϵₜₙ((1+y+ix)/2)
    Corollary 6.4 pp 22
    """
function EB(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′, a, N₀, p, U = est_utility(t,σ,T)
    # Recall that (1+y+im*x)/2 = 1-s⁺(x,y)'
    # Also note that α(1-s') = α(1-s)'
    # Hence rbstar, below is the same as the real
    # part of bstar which appears in the def of B.
    rbstar = real(1-s' + t/2*α(1-s'))
    ans = 0.0
    err = 0.0
    for n in 1:N₀
        ans += bᵗₙ(t,n)*big(e)^(-log(n)*rbstar)*ϵₜₙ(t,n,1-s')
    end
    return abs(Mₜ(t,1-s'))*ans
end

""" eB(t::Real,x::Real,y::Real)

    Returns EB/|B₀|
    """
eB(t::Real,x::Real,y::Real) = EB(t,x,y)/abs(B₀(t,x,y))

function EC(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′ = T+π*t/8  # = x/2 +πt/8 assuming σ+iT = (1+i(x+iy))/2 (66) pp 22
    # Recall (1-y+im*x)/2 = s⁺(x,y)
    # and    (1+y+im*x)/2 = 1-s⁺(x,y)'
    return big(e)^(t*π^2/64)*abs(M₀(im*T′))*(ϵ̃(t,s)+ϵ̃(t,1-s'))
end

""" eC(t::Real,x::Real,y::Real)

    Returns EC/|B₀|
    """
eC(t::Real,x::Real,y::Real) = EC(t,x,y)/abs(B₀(t,x,y))

function EC₀(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′ = T+π*t/8  # = x/2 +πt/8 assuming σ+iT = (1+i(x+iy))/2 (66) pp 22
    return big(e)^(t*π^2/64)*abs(M₀(im*T′))*(1+ϵ̃(t,s)+ϵ̃(t,1-s'))
end

""" eC₀(t::Real,x::Real,y::Real)

    Returns EC₀/|B₀|
    """
eC₀(t::Real,x::Real,y::Real) = EC₀(t,x,y)/abs(B₀(t,x,y))
                                         
        
""" ϵ̃(t::Real, σ::Real, T::Real)

    Returns a multiprecision calculation of ϵ̃, Definition (57) pp 16 of the reference.
    """
function ϵ̃(t::Real, σ::Real, T::Real)
    T′ = big(T) - π*t/8
    a = √(T′/(2*π))
    return ((0.397*9^big(σ))/(a-0.125) + 5/(3*(T′-3.33)))*big(e)^(3.49/(T′-3.33))
end

function ϵ̃(t::Real, u::Number)
    return ϵ̃(t,real(u),imag(u))
end


