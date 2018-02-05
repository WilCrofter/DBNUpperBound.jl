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
    a=b=0.0
    for n in 1:n_max
        y = π*n^2
        d = y*exp(-y*x)
        a += d
        b += y*d
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

function ΦKKL(u;n_max=100)
    return 2*ϕ(exp(2*u))
end
