
export Dirichlet_convolution, mollifiers, αterm
export bound77

using Primes

function Dirichlet_convolution(f̂::Vector{T1}, ĝ::Vector{T2}) where {T1<:Number, T2<:Number}
    N=length(f̂)
    D=length(ĝ)
    ans = zeros(promote_type(T1,T2),N*D)
    for a in 1:N, b in 1:D
        ans[a*b] += f̂[a]*ĝ[b]
    end
    return ans
end

function mollifiers(P::Integer, t::Real)
    primes = find([isprime(n) for n in 1:P]) # primes ≤ P
    b = [bᵗₙ(t,p) for p in primes] # coefficients in the Euler product 
    D = prod(primes) # number of mollifier coefficients including 0's
    λ = zeros(BigFloat,D)
    k=length(primes) # log₂(number of square-free integers in [1,D])
    for i in 0:(2^k-1) # for each pattern of exponents of a square-free product of primes ≤ P
        ibase2 = base(2,i,k) # base 2 representation of i as a string padded to k characters
        idx = find([ibase2[j]=='1' for j in 1:k])  # indices of primes involved in the product
        λ[prod(primes[idx])] = prod(b[idx])*(-1)^length(idx) # mollifier coefficient
    end
    return λ
end


""" bound77(t::Real, x::Real, y::Real)

    Returns a lower bound to |fₜ(x+iy)| formally equivalent to that given by inequality (77) pp. 29.
    Inequality (77) is |fₜ(x+iy)| ≥ 1 - |γ| - ∑₂ bᵗₙ(1+|γ|n^(y-ℜ(κ)))/n^σ
    The summation terms are formally equivalent to the absolute values of terms which define fₜ,
    bᵗₙ/n^(1-s+tα(1-s))/2) and γ⋅bᵗₙ/(s+tα(s))/2).
    """
function bound77(t::Real, x::Real, y::Real)
    absγ = abs(γₜ(t,x,y))
    s=s⁺(x,y)
    bstar=1-s+t/2*α(1-s)
    astar=s+t/2*α(s)
    bound = 1-absγ
    for n in 2:N(t,x)
        bound -= abs(Bterm(t,n,x,y,s=s,bstar=bstar))+abs(Aterm(t,n,x,y,s=s,astar=astar))
    end
    return bound
end
    
                          
