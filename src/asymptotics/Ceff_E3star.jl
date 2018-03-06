
using DBNUpperBound
using DBNUpperBound.Asymptotics

""" 

    For documentation, see the Ceff_E3star notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/Ceff_E3star.ipynb.
    """
function E3star(t::Ty1, s::Ty3) where {Ty1<:Real, Ty3<:Number}
    T0 = imag(s)
    T0 ≥ 100 || error("T0 = ",T0," is too small. T0 ≥ 100 is required for bounds to hold.")
    T0′ = T0 + π*t/8
    a0 = √(T0′/(2*π))
    return abs(1/8*√(π)*bigexp(-t*π^2/64 - π*T0/4 + 0.181/(T0′-3.33))*T0′^(3/2)*(1+5.15/(a0-1.25)))
end

""" 

    For documentation, see the Ceff_E3star notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/Ceff_E3star.ipynb.
    """
function Ceff(t::Ty1, s::Ty2) where {Ty1<:Real, Ty2<:Number}
    T′ = imag(s)+π*t/8
    s′ = s + im*π*t/8
    a  = √(T′/(2*π))
    N  = floor(Int,a)
    p  = 1 - 2*(a-N)
    U  = exp(-im*(T′/2*log(T′/(2*π)) - T′/2 - π/8))
    C0 = exp(π*im*(p^2/2+3/8) - im*√(2)*cos(π*p/2))/(2*cos(π*p))
    σ  = real(s′)
    # This precaution will surely be unnecessary in practice but...
    if !(s′≈ 0) && !(s′ ≈ 1)
        return 1/8*exp(t*π^2/64)*s′*(s′-1)/2*(-1)^N*(π^(-s′/2)*Γ(s′/2)*a^(-σ)*C0*U + π^(-(1-s′)/2)*Γ((1-s′)/2)*a^(-(1-σ))*C0'*U')
    elseif s′≈ 0
        # s′/2 ≈ 0 is a pole of Γ, but s′/2 * Γ(s′/2) ≈ 1
        return 1/8*exp(t*π^2/64)*(s′-1)*(-1)^N*(π^(-s′/2)*1*a^(-σ)*C0*U + s′/2*π^(-(1-s′)/2)*Γ((1-s′)/2)*a^(-(1-σ))*C0'*U')
    else
        # (s′-1) ≈ 0 is a pole of Γ, but (s'-1)/2 * Γ((1-s′)/2) ≈ -1
        return 1/8*exp(t*π^2/64)*s′*(-1)^N*((s′-1)/2*π^(-s′/2)*Γ(s′/2)*a^(-σ)*C0*U + π^(-(1-s′)/2)*(-1)*a^(-(1-σ))*C0'*U')
    end
end
