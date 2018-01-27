""" Φ_term(u, n; PI=π)

    Returns the value of the nth term of the infinite series, Φ(u). (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).

    """
function Φ_term{T<:Real}(u::T, n::Int; PI=convert(promote_type(T,Float64),π))
    x = PI*n^2*exp(4*u)
    return (2*PI^2*n^4*exp(9*u-x)-3*PI*n^2*exp(5*u-x))
end

""" Φ(u; n_max=100)
    
    Returns a partial summation (default 100 terms) of the infinite series, Φ(u). (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).

    A value of π may be provided using keyword, PI. The default value is Julia's value of π computed to the same precision as that of u. Typically u would be Float64 or BigFloat.

    """
function Φ{T<:Real}(u::T; n_max::Int=100, PI=convert(promote_type(T,Float64),π))
    running_sum = T(0)
    for n in 1:n_max
        running_sum += Φ_term(u,n,PI=PI)
    end
    return running_sum
end

""" phi_decay
    Alias for Φ, to match KM's Python notation
    """
phi_decay = Φ #alias

""" KM_PI
    20-digit approximation to π used in KM's Python code.
    """
KM_PI = 3.14159265358979323846
