
export Dirichlet_convolution, mollifiers, αβ
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

""" αβ(t::Real, x::Real, y::Real)

    returns vectors of coefficients, αₙ = γ*bᵗₙ*n^(y-κ̄) and βₙ=bᵗₙ, as per expression above (76) pp. 27.
    """
function αβ(t::Real, x::Real, y::Real)
    k = κ(t,x,y)
    absγ = abs(γₜ(t,x,y))
    return [absγ*bᵗₙ(t,n)*big(e)^(log(n)*(y-k')) for n in 1:N(t,x)], [bᵗₙ(t,n) for n in 1:N(t,x)]
end




""" bound77(t::Real, x::Real, y::Real)

    Returns a lower bound to |fₜ(x+iy)| formally equivalent to that given by inequality (77) pp. 29

    
    """
function bound77(t::Real, x::Real, y::Real)
    σ = real(sstar(t,x,y))
    k = real(κ(t,x,y))
    absγ = abs(γₜ(t,x,y))
    bound = 1-absγ
    for n in 2:N(t,x)
        bound -= bᵗₙ(t,n)*(1+absγ*big(e)^(log(n)*(y-k)))*big(e)^(-log(n)*σ)
    end
    return bound
end
    
                          
