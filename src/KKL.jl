# ϕ ψ π Φ

function ψ(x; n_max=100)
    ans = 0.0
    for n in 1:n_max
        ans += exp(-n^2*π*x)
    end
    return(ans)
end

function ψx(x; n_max=100)
    ans = 0.0
    for n in 1:n_max
        ans += -exp(-n^2*π*x+2*log(n)+log(π))
    end
    return(ans)
end

function ψxx(x; n_max=100)
    ans = 0.0
    for n in 1:n_max
        ans += exp(-n^2*π*x+4*log(n)+2*log(π))
    end
    return(ans)
end

function ϕ4test(x;n_max=100)
    a=3*ψx(x)
    b=2*x*ψxx(x)
    if -a > b
        return x^(1.25)*b*(1+a/b)
    elseif b > -a
        return x^(1.25)*a*(b/a+1)
    else 
        return x^(1.25)*(a+b)
    end
end

function ϕ(x; n_max=100)
    a=b=d=0.0
    for n in 1:n_max
        y = π*n^2
        d = y*exp(-y*x)
        a += d
        b += y*d
    end
    if d/a > 1e-6 || d/b > 1e-6
        warn("The partial sum has not converged. n_max should be increased")
    end
    a *= -3
    b *= 2*x
    if -a > b
        ans = x^(1.25)*b*(1+a/b)
    elseif b > -a
        ans = x^(1.25)*a*(b/a+1)
    else 
        ans = x^(1.25)*(a+b)
    end
    return ans
end

function Ξ_integrand(λ,u,z;n_max=100)
    lusq = λ*u^2
    phi = Φ(u,n_max=n_max)
    if -Inf < lusq < Inf && phi ≈ 0.0
        return 0.0
    else
        return exp(lusq)*phi*cos(u*z)
    end
end

function Ξ(λ,z;n_max=100, upper_limit = 10.0)
    ans, err = quadgk((u)-> Ξ_integrand(λ,u,z;n_max=n_max), 0.0, upper_limit, abstol=eps(λ), maxevals=10^7)
    return ans/2, err
end

function Φ(u;n_max=100)
    return ϕ(exp(2*u),n_max=n_max)
end
