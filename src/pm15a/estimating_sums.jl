
export NᵢRange, FN, FN_tail, EBound

using SpecialFunctions

#=
The following is based on http://michaelnielsen.org/polymath1/index.php?title=Estimating_a_sum, specifically the subsection "Estimating many sums."

Of interest are sums of the form F(N,z,w) = ∑ᴺ n^(-z+w⋅log(n)) where z=z₀+ζ and where ζ and w are O(1). The sum is to be estimated in terms of an exact initial sum, F(N₀,z,w), where N₀<N and a estimated partial sums over intervals of length H+1 thereafter.

A parameter, T, the length of a partial Taylor series, is required in addition to parameters already mentioned.


=#

""" hRange(H::Integer)

    Return a range of indices of length ≈H convenient for computing estimations given in http://michaelnielsen.org/polymath1/index.php?title=Estimating_a_sum.
    """
function hRange(H::Integer)
    H2=ceil(typeof(H),H/2)
    return (-H2):1:(H2)
end

""" NᵢRange(N::Integer, N₀::Integer, H::Integer)

    Given N, N₀, H, return a UnitRange for the Nᵢ such that the intervals Nᵢ+h are disjoint as h varies over this range.
    """
function NᵢRange(N::Integer, N₀::Integer, H::Integer)
    0 < N₀ < N || error("We must have 0 < Nᵢ < N")
    H2 = hRange(H)[end]
    N₁ = N₀+H2+1
    H′=2*H2+1
    return N₁:H′:(N₀+H′*ceil(typeof(N),(N-N₀)/H′))
end

""" Aᵢ(ζ::Number, w::Number, Nᵢ::Integer)
    
    Return Aᵢ(ζ,w) :=  -ζ⋅log(Nᵢ) + w⋅log²(Nᵢ)
    """
function Aᵢ(ζ::Number, w::Number, Nᵢ::Integer)
    return ζ*log(Nᵢ)+w*log(Nᵢ)^2
end

""" Bᵢ(ζ::Number, w::Number, Nᵢ::Integer)
    
    Return Bᵢ(ζ,w) := -ζ + 2⋅w⋅log(Nᵢ).
    """
function Bᵢ(ζ::Number, w::Number, Nᵢ::Integer)
    return -ζ+2*w*log(Nᵢ)
end

""" δ(ζ::Number, w::Number, N::Integer, N₀::Integer, H::Integer)
    Return δ := (|ζ|+2⋅|w|⋅log(N)⋅log(1+H/(2N₀))+|w|⋅log²(1+H/(2N₀)) for use in bounding the error term.
    """
function δ(ζ::Number, w::Number, N::Integer, N₀::Integer, H::Integer)
    leps = log(1.0+H/(2*N₀))
    return abs(ζ)+2*abs(w)*log(N)*leps+abs(w)*leps^2
end

""" Ebound(z₀::Number, ζ::Number, w::Number, N::Integer, N₀::Integer, H::Integer, T::Integer; 
           centers::StepRange{Int,Int}=NᵢRange(N,N₀,H))

    Return the right hand side of the inequality E ≤ δ^(T+1)/(T+1)! e^δ ∑ᵢe^(ℜ(Aᵢ(ζ,w))∑ₕ(Nᵢ+h)^(-ℜ(z₀)) where the summation over i is determined by the optional keyword argument, centers, which can be used, for instance, to compute the contribution of a subset of intervals. E.g., setting centers = NᵢRange(N,N₀,H)[4:4] will compute the contribution of the 4th interval, centers = NᵢRange(N,N₀,H)[4:6] will compute summed contributions of intervals 4, 5 and 6.
    """
function Ebound(z₀::Number, ζ::Number, w::Number, N::Integer, N₀::Integer, H::Integer, T::Integer;
                centers::StepRange{Int,Int}=NᵢRange(N,N₀,H))
    hidx=hRange(H)
    isum = big(0.0)
    bige = big(e)
    for Nᵢ in centers
        hsum = big(0.0)
        for h in hidx
            hsum += bige^(-real(z₀)*log(Nᵢ+h))
        end
        isum += bige^(real(Aᵢ(ζ,w,Nᵢ)))*hsum
    end
    d = δ(ζ,w,N,N₀,length(hidx))
    # Note, using e^(-logΓ(T+1+1)) for 1/(T+1)!
    return bige^(d+(T+1)*log(d)-lgamma(T+2))*isum
end

""" σᵢⱼ(ζ::Number, w::Number, Nᵢ::Integer, j::Integer, T::Integer)

    Return σᵢⱼ(ζ,w) := ∑ Bᵢ(ζ,w)^(j₁)⋅w^(j₂)/(j₁!⋅j₂!)
    where summation is over j₁≥0, j₂≥0 such that j₁+2j₂ = j, j₁+j₂ ≤ T
    """
function σᵢⱼ(ζ::Number, w::Number, Nᵢ::Integer, j::Integer, T::Integer)
    logBᵢ=log(Bᵢ(ζ,w,Nᵢ))
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

""" εᵢₕ(Nᵢ::Integer, h::Integer)
    Returns εᵢₕ := log(1 + h/Nᵢ)
    """
function εᵢₕ(Nᵢ::Integer, h::Integer)
    return log(1.0 + big(h/Nᵢ))
end

""" βᵢⱼₕ(z₀::Number, Nᵢ::Integer, j::Integer, H::Integer)

    βᵢⱼₕ := ∑ₕ (Nᵢ + h)^(-z₀)⋅(εᵢₕ)ʲ where h ranges between -H/2 and H/2
    """
function βᵢⱼₕ(z₀::Number, Nᵢ::Integer, j::Integer, H::Integer)
    ans = big(0)
    bige=big(e)
    for h in hRange(H)
        ans += bige^(-z₀*log(Nᵢ+h))*εᵢₕ(Nᵢ,h)^j
    end
    return ans
end

""" FN(z₀::Number, ζ::Number, w::Number, N₀::Integer)

    Returns the exact sum  F(N₀,z₀+ζ,w) = ∑ n^(-z₀-ζ+w⋅log(n)) where summation is from 1 to N₀.
    """
function FN(z₀::Number, ζ::Number, w::Number, N₀::Integer)
    head = big(0.0)
    z = z₀+ζ
    bige = big(e)
    for n in 1:N₀
        logn = log(n)
        head += bige^(-z*logn+w*logn^2)
    end
    return head
end

""" FN_tail(z₀::Number, ζ::Number, w::Number, N::Integer, N₀::Integer, H::Integer, T::Integer;
            centers = NᵢRange(N,N₀,H))

    Following the procedure given in http://michaelnielsen.org/polymath1/index.php?title=Estimating_a_sum, specifically the subsection "Estimating many sums," a sum F(N,z₀+ζ,w)=∑ n^(-z₀-ζ+w⋅log(n)) is estimated in terms of an exact initial sum, F(N₀,z,w), N₀<N, and estimated partial sums over intervals of length ≈H thereafter. Thus,

      FN(z₀+ζ, w) = FN₀(z₀+ζ, w) + ∑ᵢ ∑ⱼβᵢⱼₕ⋅σᵢⱼ(ζ,w) + O≤(E)

where i varies over intervals of length ≈H, j varies from 0 to 2T, and βᵢⱼₕ and σᵢⱼ are approximations. 

    This function returns an estimate of ∑ n^(-z₀-ζ+w⋅log(n)) and an associated error bound, where summation is over a range of intervals of length ≈H between N₀ and N. The summation may be over the entire range or a subset thereof, as determined by the keyword argument, centers, whose default is the entire range. It is assumed that N>N₀, z=z₀+ζ and ζ and w are O(1). Argument T is the length of a truncated Taylor series for eᵃ. 
    
    The estimated sum has the form, ∑ᵢ ∑ⱼβᵢⱼₕ⋅σᵢⱼ(ζ,w) (see functions βᵢⱼₕ, σᵢⱼ) where j varies from 0 to 2T and i varies over intervals determined by the optional keyword argument, centers. Thus
    """
function FN_tail(z₀::Number, ζ::Number, w::Number, N::Integer, N₀::Integer, H::Integer, T::Integer;
            centers = NᵢRange(N,N₀,H))
    bige = big(e)
    z=z₀+ζ
    tail = big(0.0)
    for Nᵢ in centers
        jsum = big(0.0)
        for j in 0:(2*T)
            jsum += βᵢⱼₕ(z₀,Nᵢ,j,H)*σᵢⱼ(ζ, w, Nᵢ, j, T)
        end
        tail+=jsum
    end
    return tail, Ebound(z₀, ζ, w, N, N₀, H, T; centers=centers)
end


    
