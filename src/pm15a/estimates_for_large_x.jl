
export ϵₜₙ, rₜₙ, RtN, A, B, C, C₀, EA, EB, EC, EC₀
export ϵ̃, eA, eB, eC, eC0
export ϵ̃ₜₙ, ẽA, ẽB

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
    return p ≈ 0.5 || p ≈ -0.5 ? (1-im)/4 :(exp(π*im*(p^2/2 + 3/8)) - im*√(2)*cos(π*p/2))/2*cos(π*p)
end

""" A(t::Real, x::Real, y::Real)

    A(x+iy) = Mₜ((1-y+ix)/2)*∑bᵗₙ/n^((1-y+ix)/2 +t/2*α((1-y+ix)/2))
    Corollary 6.4, pp 22. 
    """
function A(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    # Recall that (1-y+im*x)/2 = s⁺(x,y) hence T=x/2 and T′= x/2+π*t/8 as required
    T′, a, N₀, p, U = est_utility(t,σ,T)
    ans = big(0.0)
    for n in 1:N₀
        # again, note that (1-y+im*x)/2 = s⁺(x,y)
        ans += bᵗₙ(t,n)*big(e)^(log(n)*(s+t/2*α(s)))
    end
    return Mₜ(t,s)*ans
end

function A(t::Real, z::Number)
    return A(t, real(z), imag(z))
end

"""
    B(x+iy) = Mₜ((1+y-ix)/2)*∑bᵗₙ/n^((1+y-ix)/2 +t/2*α((1+y-ix)/2))
    Corollary 6.4, pp 22. 
    """
function B(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    # Recall that (1-y+ix)/2 = s⁺(x,y) hence T=x/2 and T′= x/2+π*t/8 as required
    T′, a, N₀, p, U = est_utility(t,σ,T)
    ans = 0.0
    # Note that (1+y-ix)/2 = 1-s⁺(x,y)
    for n in 1:N₀
        ans += bᵗₙ(t,n)*big(e)^(log(n)*(1-s+t/2*α(1-s)))
    end
    return Mₜ(t,1-s)*ans
end


function C(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′, a, N₀, p, U = est_utility(t,σ,T)
    return 2*(-1)^N₀*big(e)^(-im*π*y/8 +t*π^2/64)*real(M₀(im*T′)*C₀(p)*U*exp(π*im/8))
end

""" EA(t::Real,x::Real,y::Real)

    EA(x+iy) :=  |Mₜ((1-y+ix)/2)|*∑bᵗₙ/n^((1-y)/2 +t/2*ℜ(α((1-y+ix)/2)))*ϵₜₙ((1-y+ix)/2)
    Corollary 6.4 pp 22
    """
function EA(t::Real,x::Real,y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′, a, N₀, p, U = est_utility(t,σ,T)
    # Recall that (1-y+im*x)/2 = s⁺(x,y)
    ans = big(0.0)
    for n in 1:N₀
        ans += bᵗₙ(t,n)*big(e)^(log(n)*((1-y)/2)+t/2*real(α(s)))*ϵₜₙ(t,n,s)
    end
    return abs(Mₜ(t,s))*ans
end

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
    ans = 0.0
    for n in 1:N₀
        ans += bᵗₙ(t,n)*big(e)^(log(n)*((1+y)/2)+t/2*real(α(1-s')))*ϵₜₙ(t,n,1-s')
    end
    return abs(Mₜ(t,1-s'))*ans
end

function EC(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′, a, N₀, p, U = est_utility(t,σ,T)
    # Recall (1-y+im*x)/2 = s⁺(x,y)
    # and    (1+y+im*x)/2 = 1-s⁺(x,y)'
    return big(e)^(t*π^2/64)*abs(M₀(im*T′))*(ϵ̃(t,s)+ϵ̃(t,1-s'))
end

function EC₀(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    s = s⁺(x,y)
    σ, T = real(s), imag(s)
    T′, a, N₀, p, U = est_utility(t,σ,T)
    return big(e)^(t*π^2/64)*abs(M₀(im*T′))*(1+ϵ̃(t,s)+ϵ̃(t,1-s'))
end
                                         
        
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

#= As of commit 0a0f729 (teorth/dbn_upper_bound) the pdf used s rather than s_subscript_* in (69) and (70). This appears to be a typo: s_star() is used here.
 =#


""" eA(t::Real, x::Real, y::Real)

    eA(x+iy) := |γ|∑nʸ*bᵗₙ/n^ℜ(s+κ)*ϵₜₙ((1-y+ix)/2)
              = |γ|∑bᵗₙ*e^((y-ℜ(s+k))log(n))*ϵₜₙ((1-y+ix)/2)
    where s=(1-y+ix)/2
    Definition (69) pp. 23.
    """
function eA(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    #= Recall (1-y+im*x)/2 = s⁺(x,y) =#
    s = s⁺(x,y)
    ans = 0.0
    skᵣ = real(s+κ(t,x,y))
    for n in 1:N(t,x)
        ans += bᵗₙ(t,n)*big(e)^((y-real(s+κ(t,x,y)))*log(n))*ϵₜₙ(t,n,s)
    end
    return abs(γₜ(t,x,y))*ans
end

function eA(t::Real, z::Number)
    return eA(t, real(z), imag(z))
end

""" eB(t::Real, x::Real, y::Real)

    eB(x+iy) := ∑nʸ*bᵗₙ/n^ℜ(s)*ϵₜₙ((1+y+ix)/2)
              = ∑bᵗₙ*e^(-ℜ(s)log(n))*ϵₜₙ((1+y+ix)/2)
    where s=(1-y+ix)/2
    Definition (70) pp. 23.
    """
function eB(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    #= Recall
    (1+y+im*x)/2 = 1-s⁺(x,y)'   
    =#
    s= s⁺(x,y)
    sᵣ = real(s_star(t,1-s'))
    ans = 0.0
    for n in 1:N(t,x)
        ans += bᵗₙ(t,n)*big(e)^(-sᵣ*log(n))*ϵₜₙ(t,n,1-s')
    end
    return ans
end

function eB(t::Real, z::Number)
    return eB(t,real(z),imag(z))
end

""" eC(t::Real, x::Real, y::Real)

    eC := exp(tπ²/64)|M₀(iT′)|(ϵ̃((1-y+ix)/2)+ϵ̃((1+y+ix)/2))/|Mₜ((1+y+ix)/2)|
    Definition (71) pp. 23.
    """
function eC(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    #= Recall
    (1-y+im*x)/2 = s⁺(x,y)
    (1+y+im*x)/2 = 1-s⁺(x,y)'
    =#
    s = s⁺(x,y)
    T = imag(s)
    T′ = T + π*t/8 # Def (66), pp 22
    return big(e)^(t*π^2/64)*abs(M₀(im*T′))*(ϵ̃(t,s) + ϵ̃(t,1-s'))/abs(Mₜ(t, 1-s'))
end

function eC(t::Real, z::Number)
    return eC(t,real(z),imag(z))
end


""" eC0(t::Real, x::Real, y::Real)

    eC0(x+iy) := exp(tπ²/64)*|M₀(iT′)|*(1+ϵ̃((1-y+ix)/2)+ϵ̃((1+y+ix)/2))/|Mₜ((1+y+ix)/2)|
    Definition (72) pp. 23.
    """
function eC0(t::Real, x::Real, y::Real)
    in_region_5(t,x,y) || error("Parameters are not in region (5)")
    #= Recall
    (1-y+im*x)/2 = s⁺(x,y)
    (1+y+im*x)/2 = 1-s⁺(x,y)'
    =#
    s = s⁺(x,y)
    T = imag(s)
    T′ = T + π*t/8 # Def (66), pp 22
    return big(e)^(t*π^2/64)*abs(M₀(im*T′))*(1 + ϵ̃(t,s) + ϵ̃(t,1-s'))/abs(Mₜ(t,1-s'))
end

function eC0(t::Real, z::Number)
    return eC0(t,real(z),imag(z))
end

function ebound_util(n::Int, t::Real, x::Real, y::Real; sᵣ::Real=real(s_star(t,x,y)))
    return bᵗₙ(t,n)*big(e)^(-sᵣ*log(n))*(big(e)^((t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66))-1)
end

"""
    Upper bound for ϵₜₙ as given in the last line of the proof of Prop 6.6 iv.
    ϵₜₙ((1±y+ix)/2) ≤ exp((t²/32*log²(x/(4πn²))+0.313)/(T-3.33))-1
    equivalently, ϵₜₙ((1±y+ix)/2) ≤ exp((t²/16*log²(x/(4πn²))+0.626)/(x-6.66))-1
    since T=x/2 by convention.

    Note that x/2-6.66 rather than x-6.66 appears in the denominator 
    """
function ϵ̃ₜₙ(t::Real, n::Int, x::Real, y::Real)
    return big(e)^((t^2/16*log(x/(4*π*n^2))^2+0.626)/(x-6.66)) - 1
end

function ϵ̃ₜₙ(t::Real, n::Int, s::Number)
    return  ϵ̃ₜₙ(t,n,real(s),imag(s))
end
 

""" 
    Returns an upper bound to eA as given in Propositon 6.6 iv, pp 23
    eA ≤ |γ|N^|κ|∑(nʸbᵗₙ/n^(ℜ(s)))*(exp((t²/16*log²(x/(4πn²))+0.626)/(x-6.66))-1)
    equivalently eA ≤ |γ|N^|κ|∑(nʸbᵗₙ/n^(ℜ(s)))*ϵ̃ₜₙ
    """
function ẽA(t::Real, x::Real, y::Real)
    bound = 0.0
    sᵣ = real(s_star(t,x,y))
    N₀=N(t,x)
    for n in 1:N₀
        bound += bᵗₙ(t,n)*big(e)^(log(n)*(y-sᵣ))*ϵ̃ₜₙ(t,n,x,y) 
    end
    return abs(γₜ(t,x,y))*big(e)^(abs(κ(t,x,y)*log(N₀)))*bound
end

""" 
    Returns an upper bound to eB as given in Propositon 6.6 v, pp 23
    """
function ẽB(t::Real, x::Real, y::Real)
    bound = 0.0
    s = s⁺(x,y)
    sᵣ = real(s_star(t,1-s'))
    N₀=N(t,x)
    for n in 1:N₀
        bound +=  bᵗₙ(t,n)*big(e)^(-sᵣ*log(n))*ϵ̃ₜₙ(t,n,s)
    end
    return bound
end


