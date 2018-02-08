""" Φpm_term(u, n)

    Returns the value of the nth term of the infinite series, Φ(u), as defined on the PolyMath 15 page http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant.

    """
function Φpm_term(u::T, n::Int)::promote_type(T,Float64) where {T<:Number}
    x = π*n^2*exp(4*u)
    return (2*π^2*n^4*exp(9*u-x)-3*π*n^2*exp(5*u-x))
end

""" Φpm(u; n_max=100)
    
    Returns a partial summation (default 100 terms) of the infinite series, Φ(u), as defined on the PolyMath 15 page http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant.

    """
function Φpm(u::T; n_max::Int=100)::promote_type(T,Float64) where {T<:Number}
    -π/8 < imag(u) < π/8 || error("The imaginary part of u must be less than π/8 in absolute value.")
    running_sum = convert(promote_type(T,Float64),0.0)
    for n in 1:n_max
        running_sum += Φpm_term(u,n)
    end
    return running_sum
end

""" phi_decay
    Alias for Φpm, to match KM's Python notation
    """
phi_decay = Φpm #alias
