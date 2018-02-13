""" BigNum

    Return type for asymptotic functions.
    """
const BigNum = Union{BigFloat, Complex{BigFloat}}

function bigify(s::T)::BigNum where {T<:Number}
    return imag(s)==0 ? big(s) : big(real(s))+im*big(imag(s))
end

function Ψpm(α::T)::BigNum where {T<:Number}
    α = bigify(α)
    return 2*π*exp(im*(π/2)*α^2-5*π/8)*cos(π*(α^2/2-α-π/8))/cos(π*α)
end

function Rpm(t::T1,N::I1,M::I2,s::T2)::BigNum where {T1<:Number, I1<:Integer, I2<:Integer, T2<:NotBigComplex}
    # Note that ξ(s)= π^(-s/2)*(s/2)*Γ(s/2)*(s-1)*ζ(s)
    bigs=bigify(s)
    return -ξ(s)*exp(-t^2/64-im*π*bigs)*(2*π*im*M)^(bigs-1)*Ψpm(s/(2*π*im*M)-N)/(2*π*im)
end

function Fpm(t::T1,N::I1,s::NBC)::BigNum where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    if s ≈ 0
        Γfactor = 1.0
    else
        Γfactor = (s/2)*Γ(s/2)
    end
    bigs=bigify(s)
    psum = 0
    for n in 1:N
        psum += exp(t*log((bigs+4)/(2*π*n^2))^2/16)/n^bigs
    end
    return (s-1)*π^(-bigs/2)*Γfactor*psum
end

function ξ(t,s,N,M)
    return Fpm(t,N,s)+Fpm(t,M,(1-s)')'+Rpm(t,N,M,s)
end
