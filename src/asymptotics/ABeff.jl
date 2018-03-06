
""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
Nseff(s,t) = floor(Int, √((imag(s) - π*t/8)/(2*π)))

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
αeff(s) = 1/(2*s) + 1/(s-1) + 1/2*log(s/(2*π))

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
H01(s) = (s*(s-1)/2)*√(2*π)*bigexp(-s/2*log(π) + (s/2-1/2)*log(s/2) - s/2)

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
function aeff(s::T1, t::T2, n::Int) where {T1<:Number, T2<:Real}
    return bigexp(-log(n)*(s + t*αeff(s)/2 - t/4*log(n)))
end

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
function beff(s::T1, t::T2, n::Int) where {T1<:Number, T2<:Real}
    return bigexp(-log(n)*( 1-s + t*αeff((1-s)')'/2 - t/4*log(n)))
end

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
function Aeff_μ(s::T1, t::T2) where {T1<:Number, T2<:Real}
    return 1/8*bigexp(t/4*αeff(s)^2)*H01(s)
end

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
function Beff_μ(s::T1, t::T2) where {T1 <: Number, T2 <: Real}
    return 1/8*bigexp(t/4*(αeff((1-s)')')^2)*H01((1-s)')'
end

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
function Aeff(t::T1, s::T2, N::Int) where {T1 <: Real, T2 <:Number}
    psum = 0.0
    for n in 1:N
        psum += aeff(s,t,n)
    end
    return Aeff_μ(s,t)*psum
end

""" 

    For documentation, see the ABeff notebook https://github.com/WilCrofter/DBNUpperBound.jl/blob/master/notebooks/ABeff.ipynb.
    """
function Beff(t::T1, s::T2, N::Int) where {T1 <: Real, T2 <:Number}
    psum = 0.0
    for n in 1:N
        psum += beff(s,t,n)
    end
    return Beff_μ(s,t)*psum
end
