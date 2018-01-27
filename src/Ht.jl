""" Ht_term(t, u, n; PI=π)

    Returns the value of the nth term of the infinite series, Φ(u)*exp(t*u^2). (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).
    """
function Ht_term{T1<:Real, T2<:Real}(t::T1, u::T2, n::Int; PI=convert(promote_type(T1,T2,Float64),π))
    x = t*u^2-PI*n^2*exp(4*u)
    return (2*PI^2*n^4*exp(9*u+x)-3*PI*n^2*exp(5*u+x))
end

function Ht_series{T1<:Real, T2<:Real, T3<:Number}(t::T1, u::T2, z::T3;
                                                   n_max=100,
                                                   PI=convert(promote_type(T1,T2,typeof(real(z)),typeof(imag(z)),Float64),π))
    PI
    ans = convert(promote_type(Complex{T1},Complex{T2},T3),0)
    for n in 1:n_max
        ans += Ht_term(t,u,n,PI=PI)
    end
    return ans
end

""" Ht(t, z; n_max=100, PI=π)
    
    Returns an approximate value of Ht(z), the integral from 0 to ∞ of Φ(u)*exp(t*u^2)*cos(uz), where the first two factors are approximated by truncated series of n_max (default 100) terms. (See http://michaelnielsen.org/polymath1/index.php?title=De_Bruijn-Newman_constant).
    """
function Ht{T1<:Real,T2<:Number}(t::T1, z::T2; n_max=100,
                                 PI=convert(promote_type(T1,typeof(real(z)),typeof(imag(z)),Float64),π))
    return quadgk((u)-> Ht_series(t,u,z; n_max=n_max, PI=PI), 0.0, Inf)
end
