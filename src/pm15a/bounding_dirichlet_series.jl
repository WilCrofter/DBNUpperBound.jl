
export Dirichlet_convolution, mollifiers

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
