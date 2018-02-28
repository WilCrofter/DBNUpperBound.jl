

function Aprime_a(t::T1, s::T2, n::Int) where {T1<:Real, T2<:Number}
    return bigexp((t/16)*log((s+4)/(2*π*n^2))^2 - s*log(n))
end

function Aprime_μ(t::T1, s::T2) where {T1<:Real, T2<:Number}
    return (2/8)*π^(-s/2)*√(2*π)*bigexp(((s+4)/2 - 1/2)*log((s+4)/2)-(s+4)/2)
end

function Aprime(t::T1, s::T2, N::Int) where {T1<:Real, T2<:Number}
    asum = 0.0
    for n = 1:N
        asum += Aprime_a(t,s,n)
    end
    return Aprime_μ(t,s)*asum
end

Bprime(t,s,N) = Aprime(t,1-s,N) # (This is legal function definition syntax in Julia)
B0prime(t,s) = Bprime(t,s,1)
