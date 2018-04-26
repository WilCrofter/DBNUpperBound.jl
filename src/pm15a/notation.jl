""" complex_power(z::Number, w::Number)

    Return exp(wLog(z)) where Log refers to the standard branch of the log (which is Julia's default.)
    """
function complex_power(z::Number, w::Number)
    !(imag(z)≤0) || error("Complex power is undefined for imag(z)≤0")
    return bigexp(w*log(z))
end


""" gaussian_identity(a,b,c)

    Return ∫exp(-(au²+bu+c))du = √(π/a)exp(b²/(4a)-c) where integration is over the real line
    """
function gaussian_identity(a::Number, b::Number, c::Number)
    return complex_power(π/a,1/2)*bigexp(b^2/(4*a) - c)
end
