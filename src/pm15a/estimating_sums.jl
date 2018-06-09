
export Nᵢ, Aᵢ

using SpecialFunctions

#=
The following is based on http://michaelnielsen.org/polymath1/index.php?title=Estimating_a_sum, specifically the subsection "Estimating many sums."

Of interest are sums of the form F(N,z,w) = ∑ᴺ n^(-z+w⋅log(n)) where z=z₀+ζ and where ζ and w are O(1). The sum is to be estimated in terms of an exact initial sum, F(N₀,z,w), where N₀<N and a estimated partial sums over intervals of length H+1 thereafter.

A parameter, T, the length of a partial Taylor series, is required in addition to parameters already mentioned.


=#

""" Nᵢ(N::Integer, N₀::Integer, H::Integer)
    Given N, N₀, H, return a UnitRange for the Nᵢ
    Note that for the intervals to be disjoint we must have Nᵢ₊₁-H/2 = 1 + Nᵢ + H/2, so  Nᵢ₊₁=1+H+Nᵢ with N₁=1+N₀.
    """
function Nᵢ(N::Integer, N₀::Integer, H::Integer)
    0 < N₀ < N || error("We must have 0 < Nᵢ < N")
    H2 = ceil(typeof(H),H/2)
    N₁ = N₀+H2+1
    H′=2*H2+1
    return N₁:H′:(N₀+H′*ceil(typeof(N),(N-N₀)/H′))
end

""" Aᵢ(ζ::Real, w::Real, Nᵢ::Integer)
    
    Return Aᵢ(ζ,w) :=  -ζ⋅log(Nᵢ) + w⋅log²(Nᵢ)
    """
function Aᵢ(ζ::Real, w::Real, Nᵢ::Integer)
    return ζ*log(Nᵢ)+w*log(Nᵢ)^2
end

""" Bᵢ(ζ::Real, w::Real, Nᵢ::Integer)
    
    Return Bᵢ(ζ,w) := -ζ + 2⋅w⋅log(Nᵢ).
    """
function Bᵢ(ζ::Real, w::Real, Nᵢ::Integer)
    return -ζ+2*w*log(Nᵢ)
end

""" δ(ζ::Real, w::Real, N::Integer, N₀::Integer, H::Integer)
    Return δ := (|ζ|+2⋅|w|⋅log(N)⋅log(1+H/(2N₀))+|w|⋅log²(1+H/(2N₀)) for use in bounding the error term.
    """
function δ(ζ::Real, w::Real, N::Integer, N₀::Integer, H::Integer)
    leps = log(1.0+H/(2*N₀))
    return abs(ζ)+2*abs(w)*log(N)*leps+abs(w)*leps^2
end

""" Ebound(z₀::Number, ζ::Real, w::Real, N::Integer, N₀::Integer, H::Integer)

    Return the right hand side of the inequality E ≤ δ^(T+1)/(T+1)! e^δ ∑ᵢe^(ℜ(Aᵢ(ζ,w))∑ₕ(Nᵢ+h)^(-ℜ(z₀))
    """
function Ebound(z₀::Number, ζ::Real, w::Real, N::Integer, N₀::Integer, H::Integer)
    Nᵢs = Nᵢ(N,N₀,H) 
    isum = big(0.0)
    bige = big(e)
    H2=round(typeof(H),H/2)
    for Nᵢ′ in Nᵢs
        hsum = big(0.0)
        for h in (Nᵢ′-H2):1:(Nᵢ′+H2)
            hsum += bige^(-real(z₀)*log(Nᵢ′+h))
        end
        isum += bige^(real(Aᵢ(ζ,w,Nᵢ′)))*hsum
    end
    d = δ(ζ,w,N,N₀,H)
    # Note, using e^(-logΓ(T+1+1)) for 1/(T+1)!
    return bige^(d+(T+1)*log(T)-lgamma(T+2))*isum
end

"""

    Return σᵢⱼ(ζ,w) := ∑ Bᵢ(ζ,w)^(j₁)⋅w^(j₂)/(j₁!⋅j₂!)
    where summation is over j₁≥0, j₂≥0 such that j₁+2j₂ = j, j₁+j₂ ≤ T
    """
function σᵢⱼ(ζ::Real, w::Real, Nᵢ′::Integer, j::Integer, T::Integer)
    logBᵢ=log(Bᵢ(ζ,w,Nᵢ′))
    logw =log(w)
    ans=big(0.0)
    bige=big(e)
    # j₁=j-2j₂ and j₁+j₂ ≤ T imply j₂≥(j-T)₊
    # j₁=j-2j₂ implies j₂≤j/2. Thus,
    for j₂ in max(j-T,0):1:floor(typeof(j),j/2)
        j₁=j-2*j₂
        ans += bige^(j₁*logBᵢ+j₂*logw-lgamma(j₁+1)-lgamma(j₂+1))
    end
    return ans
end

"""
    Returns εᵢₕ := log(1 + h/N_i)
    """
function εᵢₕ(Nᵢ′::Integer, h::Integer)
    return log(1.0 + h/Nᵢ′)
end

""" 

    βᵢⱼₕ := ∑ₕ (Nᵢ + h)^(-z₀)⋅(εᵢₕ)ʲ where h ranges between -H/2 and H/2
    """
function βᵢⱼₕ(z₀::Number, Nᵢ′::Integer, H::Integer)
    ans = big(0)
    bige=big(e)
    H2 = ceil(typeof(H),H/2)
    for h in -H2:1:(H-H2)
        ans += bige^(-z₀*log(Nᵢ′+h)+j*log(εᵢⱼ(Nᵢ′+h)))
    end
    return ans
end
