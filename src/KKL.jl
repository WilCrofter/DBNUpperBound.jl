# ϕ ψ π Φ

function ψKKL(x::T; n_max=100) where {T<:Number}
    real(x) > 0 || error("For convergence, x must have positive real part")
    ans = 0.0
    for n in 1:n_max
        ans += exp(-n^2*π*x)
    end
    return(ans)
end


function ψxKKL(x::T; n_max=100) where {T<:Number}
    real(x) > 0 || error("For convergence, x must have positive real part")
    ans = 0.0
    for n in 1:n_max
        ans += -n^2*exp(-n^2*π*x)
    end
    return(ans*π)
end

function ψxxKKL(x::T; n_max=100) where {T<:Number}
    real(x) > 0 || error("For convergence, x must have positive real part")
    ans = 0.0
    for n in 1:n_max
        ans += n^4*exp(-n^2*π*x)
    end
    return(ans*π^2)
end

function ϕKKL4test(x::T; n_max=100) where {T<:Number}
    return x^(5/4)*(2*x*ψxxKKL(x; n_max=n_max) + 3*ψxKKL(x; n_max=n_max))
end

function ϕKKL(x::T; n_max=100) where {T<:Number}
    real(x) > 0 || error("For convergence, x must have positive real part")
    a=b=d=0.0
    for n in 1:n_max
        y = π*n^2
        d = y*exp(-y*x)
        a += d
        b += y*d
    end
    if abs(d/a) > 1e-6 || abs(d/b) > 1e-6
        warn("The partial sum in ϕ($x) has not converged. n_max should be increased")
    end
    a *= -3
    b *= 2*x
    ans = x^(1.25)*(a+b)
    return ans
end

function ΦKKL(u::T;n_max=100) where {T<:Number}
    return 2*ϕKKL(exp(2*u),n_max=n_max)
end

function ΞKKL_integrand(λ::T1,u::T2,z::T3;n_max=100) where {T1<:Real, T2<:Number, T3<:Number}
    lusq = λ*u^2
    phi = ΦKKL(u,n_max=n_max)
    if -Inf < real(lusq) < Inf && phi ≈ 0.0
        return 0.0
    else
        return exp(lusq)*phi*cos(u*z)
    end
end

function ΞKKL(λ::T1,z::T2;n_max=100, upper_limit = 10.0) where {T1<:Real, T2<:Number}
    ans, err = quadgk((u)-> ΞKKL_integrand(λ,u,z;n_max=n_max), 0.0, upper_limit, abstol=eps(λ), maxevals=10^7)
    return ans/2, err
end

