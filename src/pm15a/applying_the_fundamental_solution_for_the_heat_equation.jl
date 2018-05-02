""" Hₜ(t::Real, x::real, y::Real)

    Estimates Hₜ via equation (34), pp 10 of the reference, 1/8∫ξ((1+iz)/2+v√(π))exp(-v²)dv/√(π).

    Because Julia's implementation of Riemann's ζ does not accept multiprecision complex numbers, this function will throw errors when x exceeds 900 (approximately.) 
    
    """
function Hₜ(t::Real, x::Real, y::Real; lower_limit::Real=-10.0, upper_limit::Real=+10.0) 
    s=s⁺(x,y)
    srt = √(t)
    q(v) = ξ(s+srt*v)*bigexp(-v^2)/8/√(π)
    return quadgk(q,lower_limit,0.0,upper_limit)
end
