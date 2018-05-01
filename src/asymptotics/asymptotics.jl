module Asymptotics

export A,B,C,Nz,Ns,bigexp,bigcos,bigΓ,Ψpm
export Aprime, Bprime, B0prime
export Aeff, Beff, Ceff
export E1, E2, E3, E3star
import DBNUpperBound.NotBigComplex, DBNUpperBound.logΓ
export θ_default, Itθ, XisOK, Itθ_tail, minimum_n, series_tail_Itθ9, series_tail_Itθ5, Ht_tail
export Jtθ, Jtθ_tail, series_tail_Jtθ9, series_tail_Jtθ5, H′t_tail, logmag, ω
export ω_pds, pwquad, pwItθ

include("ABprime.jl")
include("ABeff.jl")
include("Ceff_E3star.jl")
include("EffectiveBounds.jl")
include("uniform_bounds.jl")
include("I_and_J.jl")


# Notation:
# Ht's argument is denoted z=x+iy where x,y are real.
# Arguments to A, B, and C include s=(1+iz)/2
# The case of interest is t=y=0.4, and x>>0.0.
# Specifically, x > 2e+5 is of interest.


function bigcos(s::T) where {T<:Number}
    return (bigexp(-im*s)+bigexp(im*s))/2
end

function bigΓ(s::T) where {T<:NotBigComplex}
    return bigexp(logΓ(s))
end

"""
    https://terrytao.wordpress.com/2018/02/02/polymath15-second-thread-generalising-the-riemann-siegel-approximate-functional-equation/#comment-492182
"""
function Ψpm(α::T) where {T<:Number}
    return 2*π*bigexp(im*(π/2)*(α^2-5/4))*bigcos(π*(α^2/2-α-1/8))/bigcos(π*α)
end

function C(t::T1,N::I1,M::I2,s::T2) where {T1<:Number, I1<:Integer, I2<:Integer, T2<:NotBigComplex}
    p1 = s/2 ≈ 0.0 ? 1.0 : (s/2)*bigΓ(s/2)
    p2 = 1-s ≈ 0.0 ? -1.0 : (s-1)*bigΓ(1-s)
    return -p1*p2*bigexp(-(t^2)/64 - s*(im*π+log(π)/2) + (s-1)*log(2*π*im*M)) * Ψpm(s/(2*π*im*M)-N)/(2*π*im)
end

function A(t::T1,N::I1,s::NBC) where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    Γfactor =  s ≈ 0 ? 1.0 : (s/2)*bigΓ(s/2)
    psum = 0.0
    for n in 1:N
        psum += bigexp( (t/16)*log((s+4)/(2*π*n^2))^2-s*log(n))
    end
    return (1/8)*(s-1)*Γfactor*bigexp(-log(π)*s/2)*psum
end

function B(t::T1,M::I1,s::NBC) where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    Γfactor =  (s-1)/2 ≈ 0 ? -1.0 : ((s-1)/2)*bigΓ((1-s)/2)
    psum = 0
    for m in 1:M
        psum += bigexp( (t/16)*log((5-s)/(2*π*m^2))^2 )/m^(1-s)
    end
    return (1/8)*s*Γfactor*bigexp(log(π)*(s-1)/2)*psum
end

function Nz(z)
    return floor(Int, sqrt(real(z)/(4*π)))
end

function Ns(s)
    return floor(Int, sqrt(imag(s)/(2*π)))
end


end
