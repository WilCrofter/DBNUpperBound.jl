""" Φ_term(u, n)

    Returns the value of the nth term of the infinite series, Φ(u). (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).

    """
function Φ_term{T<:Number}(u::T, n::Int)::promote_type(T,Float64)
    x = π*n^2*exp(4*u)
    return (2*π^2*n^4*exp(9*u-x)-3*π*n^2*exp(5*u-x))
end

""" Φ(u; n_max=100)
    
    Returns a partial summation (default 100 terms) of the infinite series, Φ(u). (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).

    """
function Φ{T<:Number}(u::T; n_max::Int=100)::promote_type(T,Float64)
    -π/8 < imag(u) < π/8 || error("The imaginary part of u must be less than π/8 in absolute value.")
    running_sum = convert(promote_type(T,Float64),0.0)
    for n in 1:n_max
        running_sum += Φ_term(u,n)
    end
    return running_sum
end

""" phi_decay
    Alias for Φ, to match KM's Python notation
    """
phi_decay = Φ #alias
