""" Hₜ(t::Real, x::Real, y::Real)

    Estimates Hₜ via equation (34), pp 10 of the reference, 1/8∫ξ((1+iz)/2+v√(π))exp(-v²)dv/√(π).

    Julia's implementation of Riemann's ζ has no complex multiprecision capability. Since our implementation of ξ depends on Julia's ζ, arguments to Hₜ are limited to single precision and will throw errors when x exceeds 900 (approximately.) 
    
    """
function Hₜ(t::Real, x::Real, y::Real; lower_limit::Real=-10.0, upper_limit::Real=+10.0) 
    s=s⁺(x,y)
    srt = √(t)
    q(v) = ξ(s+srt*v)*big(e)^(-v^2)/8/√(π)
    return quadgk(q,lower_limit,0.0,upper_limit)
end

function Hₜ(t::Real, z::Number; lower_limit::Real=-10.0, upper_limit::Real=+10.0)
    return Hₜ(t,real(z),lower_limit=lower_limit, upper_limit=upper_limit)
end
