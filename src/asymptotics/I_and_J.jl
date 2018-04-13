
""" θ_default(y)

    Returns a value for θ which should cancel most oscillation in Itθ or Jtθ.
    """
θ_default(y) = π/8 - 1/4*atan((abs(y)+9)/abs(y)) # function definition


""" logmag(x,y,t,σ,β,θ,m)

    Returns the log magnitude of Itθ's integrand.
    """
function logmag(x,y,t,σ,β,θ,m)
    return -e^(4*σ)*β*cos(4*θ)+m*σ-t*θ^2+t*σ^2-x*θ-y*σ
end

""" ω(x,y,t,σ,β,θ,m)

    Returns the angle of Itθ's integrand.
    """
function ω(x,y,t,σ,β,θ,m)
    return -e^(4*σ)*β*sin(4*θ)+m*θ+2*t*θ*σ+x*σ-y*θ
end

""" ω_pds(x,t,β,θ,lower_limit, upper_limit; scale=1.0)

    Returns a vector such that ω(x,y,t,σ[i+1],β,θ,m)-ω(x,y,t,σ[i],β,θ,m) ≈ 2π.
    """
function ω_pds(x,t,β,θ,lower_limit, upper_limit; scale=1.0)
    ans = [lower_limit]
    while ans[end] < upper_limit
        push!(ans, ans[end]+scale*2*π/abs(2*t*θ+x-4*e^(4*ans[end])*β*sin(4*θ)))
    end
    return ans
end

""" pwquad(f,limits; reltol=sqrt(eps(typeof(limits[1]))), abstol=0, maxevals=10^7, order=7, norm=vecnorm)

    Integrates f by piecewise Gauss-Kronrod quadrature over intervals limits[i] ≤ x ≤ limits[i+1], returning accumulated results and error estimates. 
    """
function pwquad(f,limits)
    val =0.0
    err =0.0
    for i in 2:length(limits)
        tmp = quadgk(f,limits[i-1],limits[i])
        val += tmp[1]
        err += tmp[2]
    end
    return val, err
end

""" pwItθ(x, y, n, m, limits; t=.4, θ=θ_default(y)

    Evaluates the integral defining Itθ by piecewise Gauss-Kronrod quadrature over intervals given in limits, returning a value, an error estimate, and a formal bound for the tail.
    """
function pwItθ(x,y,n,m,limits; t=.4, θ=θ_default(y))
    β = π*n^2
    f(σ) = exp(logmag(x,y,t,σ,β,θ,m))*(cos(ω(x,y,t,σ,β,θ,m))+im*sin(ω(x,y,t,σ,β,θ,m)))
    val, err = pwquad(f,limits)
    return val, err, Itθ_tail(limits[end], x, y, n, m; t=t, θ = θ)
end

""" Itθ(x, y, n, m; t=.4, lower_limit=big(0.0), upper_limit=big(Inf), θ = θ_default(y))

    Estimates Itθ by adaptive Gauss-Kronrod quadrature, returning Itθ and an error estimate.
    NOTE: GK quadrature will fail for oscillatory arguments--larger values of x.

    The quadrature function uses multiprecision arithmetic if at least one limit of integration is multiprecision.

    If a very unhelpful diagnostic appears, increase precision. It is 256 bits by default. 
    """
function Itθ(x, y, n, m; t=.4, lower_limit=big(0.0), upper_limit=big(Inf), θ = θ_default(y))
        β = π*n^2
        mag(σ) = exp(logmag(x,y,t,σ,β,θ,m))
        ω1(σ) = ω(x,y,t,σ,β,θ,m)
        return quadgk((σ)->mag(σ)*(cos(ω1(σ))+im*sin(ω1(σ))), lower_limit, upper_limit) # return quadrature, error estimate
end

""" XisOK(X, y, n, m; t=.4, θ=θ_default(y))

    Returns true iff β*bigexp(4*X)*cos(4*θ) > max(t/2, (a+2*t*X)/4) indicating bounds for the tails of Itθ and Jtθ will be valid.
    """
function XisOK(X, y, n, m; t=.4, θ=θ_default(y))
    β = π*n^2
    a = m+y
    return β*bigexp(4*X)*cos(4*θ) > max(t/2, (a+2*t*X)/4)
end

""" Itθ_tail(X, x, y, n, m; t=.4, θ = θ_default(y))
    
    Returns a provably correct upper bound for the tail of Itθ, meaning the integral from X to ∞.
    """
function Itθ_tail(X, x, y, n, m; t=.4, θ = θ_default(y))
    XisOK(X,y,n,m,t=t,θ=θ) || error("X is too small")
    β = π*n^2
    a = m+y
   return bigexp(-t*θ^2+t*X^2-β*e^(4*X)*cos(4*θ)-θ*x+a*X)/(4*β*e^(4*X)*cos(4*θ)-a-2*t*X)
end

""" minimum_n(x,y,m; t= .4, θ=θ_default(y))

    Return the minimum n for which π*n^2cos(4θ)>max(t/2,a/4), where a=m+y. Note that cos(4θ)≥0.0 for allowable θ. 
    """
function minimum_n(x,y,m; t= .4, θ=θ_default(y))
    return ceil(Int, √(max(t/2, (m+y)/4)/(π*cos(4*θ))))
end

""" series_tail_Itθ9(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,9, t=t, θ=θ))

    Return a provably correct upper bound on the tail |∑2*π^2*n^4*Itθ(x,y,9)|, along with n0 where summation is from n0 to ∞.
    """
function series_tail_Itθ9(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,9,t=t,θ=θ))
    α = bigexp(-π*n0*cos(4*θ))
    m = 9 # by definition of this series
    a = m+y
    return (2*π^2*bigexp(-t*θ^2-θ*x-π*n0^2*cos(4*θ)))/(4*π*n0^2*cos(4*θ)-a)*
            (n0^4/(1-α) + 4*n0^3*α/(1-α)^2 + 6*n0^2*(α^2 + α)/(1-α)^3 + 4*n0*(α^3+4*α^2+α)/(1-α)^4 + (α^4+11*α^3+11*α^2+α)/(1-α)^5), n0
end

""" series_tail_Itθ5(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,5, t=t, θ=θ))

    Return a provably correct upper bound on the tail |∑3*π*n^2*Itθ(x,y,5)| where summation is from n0 to ∞.
    """
function series_tail_Itθ5(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,5,t=t,θ=θ))
    α = bigexp(-π*n0*cos(4*θ))
    m = 5 # by definition of this series
    a = m+y
    return (3*π*bigexp(-t*θ^2-θ*x-π*n0^2*cos(4*θ)))/(4*π*n0^2*cos(4*θ)-a)*
           (n0^2/(1-α) + 2*n0*α/(1-α)^2 + (α^2 + α)/(1-α)^3), n0
end

""" Ht_tail(x,y; t=.4, θ=θ_default(y), n0=max(minimum_n(x,y,5,t=t,θ=θ), minimum_n(x,y,9,t=t,θ=θ))
    
    Return an upper bound for the tail of a series expressing Ht in terms of Itθ along with the index of the tail's first term.
    """
function Ht_tail(x,y; t=.4, θ=θ_default(y), n0=max(minimum_n(x,y,5,t=t,θ=θ), minimum_n(x,y,9,t=t,θ=θ)))
    return series_tail_Itθ9(x,y,t=t,θ=θ,n0=n0)[1]+series_tail_Itθ5(x,y,t=t,θ=θ,n0=n0)[1] +
               series_tail_Itθ9(x,-y,t=t,θ=θ,n0=n0)[1]+series_tail_Itθ5(x,-y,t=t,θ=θ,n0=n0)[1], n0
end

""" Jtθ(x, y, n, m; t=.4, lower_limit=big(0.0), upper_limit=big(Inf), θ = θ_default(y))

    Estimates Jtθ by adaptive Gauss-Kronrod quadrature, returning Jtθ and an error estimate. 
    NOTE: GK quadrature will fail for oscillatory arguments--larger values of x.
    
    The quadrature function uses multiprecision arithmetic if at least one limit of integration is multiprecision.

    If a very unhelpful diagnostic appears, increase precision. It is 256 bits by default. 
    """
function Jtθ(x, y, n, m; t=.4, lower_limit=big(0.0), upper_limit=big(Inf), θ = θ_default(y))
    β = π*n^2
    mag(σ) = exp(logmag(x,y,t,σ,β,θ,m))
    ω1(σ) = ω(x,y,t,σ,β,θ,m)
    return quadgk((σ)->mag(σ)*(σ*cos(ω1(σ))-θ*sin(ω1(σ))+im*(σ*sin(ω1(σ))+θ*cos(ω1(σ)))), lower_limit, upper_limit) # return quadrature, error estimate
end

""" Jtθ_tail(X, x, y, n, m; t=.4, θ = θ_default(y))
    
    Returns a provably correct upper bound for the tail of Jtθ, meaning the integral from X to ∞.
    """
function Jtθ_tail(X, x, y, n, m; t=.4, θ = θ_default(y))
    XisOK(X,y,n,m,t=t,θ=θ) || error("X is too small")
    β = π*n^2
    a = m+y
   return bigexp(-t*θ^2+t*X^2-β*e^(4*X)*cos(4*θ)-θ*x+a*X)*(abs(X+im*θ)/(4*β*e^(4*X)*cos(4*θ)-a-2*t*X) 
                                                           + 1/(4*β*e^(4*X)*cos(4*θ)-a-2*t*X)^2)
end

""" series_tail_Jtθ9(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,9, t=t, θ=θ))

    Returns a provably correct upper bound on the tail |∑2*π^2*n^4*Jtθ(x,y,9)| where summation is from n0 to ∞.
    """
function series_tail_Jtθ9(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,9,t=t,θ=θ))
    α = bigexp(-π*n0*cos(4*θ))
    m = 9 # by definition of this series
    a = m+y
    return (2*π^2*bigexp(-t*θ^2-π*n0^2*cos(4*θ)-θ*x))*(θ/(4*π*n0^2*cos(4*θ)-a) + 1/(4*π*n0^2*cos(4*θ)-a)^2)
end


""" series_tail_Jtθ5(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,5, t=t, θ=θ))

    Return a provably correct upper bound on the tail |∑3*π*n^2*Jtθ(x,y,5)| where summation is from n0 to ∞.
    """
function series_tail_Jtθ5(x,y; t=.4, θ=θ_default(y), n0=minimum_n(x,y,5,t=t,θ=θ))
    α = bigexp(-π*n0*cos(4*θ))
    m = 5 # by definition of this series
    a = m+y
    return 3*π*bigexp(-t*θ^2-π*n0^2*cos(4*θ)-θ*x)*(θ/(4*π*n0^2*cos(4*θ)-a) + 1/(4*π*n0^2*cos(4*θ)-a)^2)
end

""" H′t_tail(x,y; t=.4, θ=θ_default(y), n0=max(minimum_n(x,y,5,t=t,θ=θ), minimum_n(x,y,9,t=t,θ=θ))
    
    Return an upper bound for the tail of a series expressing H′t in terms of Jtθ along with the index of the tail's first term.
    """
function H′t_tail(x,y; t=.4, θ=θ_default(y), n0=max(minimum_n(x,y,5,t=t,θ=θ), minimum_n(x,y,9,t=t,θ=θ)))
    return series_tail_Jtθ9(x,y,t=t,θ=θ,n0=n0)[1]+series_tail_Jtθ5(x,y,t=t,θ=θ,n0=n0)[1] +
               series_tail_Jtθ9(x,-y,t=t,θ=θ,n0=n0)[1]+series_tail_Jtθ5(x,-y,t=t,θ=θ,n0=n0)[1], n0
end

