module Asymptotics

export A,B,C,Nz,Ns,bigexp,bigcos,bigΓ,Ψpm
import DBNUpperBound.NotBigComplex, DBNUpperBound.logΓ

include("ABprime.jl")

function bigify(s::T) where {T<:Number}
    return imag(s)==0 ? big(real(s)) : big(real(s))+im*big(imag(s))
end

# Notation:
# Ht's argument is denoted z=x+iy where x,y are real.
# Arguments to A, B, and C include s=(1+iz)/2
# The case of interest is t=y=0.4, and x>>0.0.
# Specifically, x > 2e+5 is of interest.

function bigexp(s::T) where {T<:Number}
    ans = exp(real(s))
    return exp(im*imag(s))*(isfinite(ans) && !(ans == 0) ? ans : exp(bigify(real(s))))
end

function bigcos(s::T) where {T<:Number}
    u,v = real(s),imag(s)
    if v < 0.0
        u=-u
        v=-v
    end
    return bigexp(v)*((1+bigexp(-2*v))*cos(u) + im*(1-bigexp(-2*v))*sin(u))
end

function bigΓ(s::T) where {T<:NotBigComplex}
    return bigexp(logΓ(s))
end

function Ψpm(α::T) where {T<:Number}
    return 2*π*bigexp(im*(π/2)*α^2-5*π/8)*bigcos(π*(α^2/2-α-π/8))/bigcos(π*α)
end

function C(t::T1,N::I1,M::I2,s::T2) where {T1<:Number, I1<:Integer, I2<:Integer, T2<:NotBigComplex}
    p1 = s/2 ≈ 0.0 ? 1.0 : (s/2)*bigΓ(s/2)
    p2 = 1-s ≈ 0.0 ? -1.0 : (s-1)*bigΓ(1-s)
    return -p1*p2*bigexp(-(t^2)/64 - s*(im*π+log(π)/2) + (s-1)*log(2*π*im*M)) * Ψpm(s/(2*π*im*M)-N)/(2*π*im)
end

function A(t::T1,N::I1,s::NBC) where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    Γfactor =  s ≈ 0 ? 1.0 : (s/2)*bigΓ(s/2)
    psum = 0
    for n in 1:N
        psum += bigexp( (t/16)*log((s+4)/(π*n^2))^2-s*log(n))
    end
    return (1/8)*(s-1)*Γfactor*bigexp(-log(π)*s/2)*psum
end

function B(t::T1,M::I1,s::NBC) where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    Γfactor =  (s-1)/2 ≈ 0 ? -1.0 : ((s-1)/2)*bigΓ((1-s)/2)
    psum = 0
    for m in 1:M
        psum += bigexp( (t/16)*log((5-s)/(π*m^2))^2 -(1-s)*log(m) )
    end
    return (1/8)*s*Γfactor*bigexp(-log(π)*s/2)*psum
end

function Nz(z)
    return floor(Int, sqrt(real(z)/(4*π)))
end

function Ns(s)
    return floor(Int, sqrt(imag(s)/(2*π)))
end


end
