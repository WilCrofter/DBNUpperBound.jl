
export Dirichlet_convolution, mollifiers, αterm
export bound77, bound78

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

""" αterm(t::Real, n::Int, x::Real, y::Real)

    Returns  αₙ/|γ| = bᵗₙ⋅n^(y-κ̄) as per expression above (76) pp. 27,
    and a term of magnitude 1, such that the product of returned values
    formally equals bᵗₙ/n^(s+t/2*α(s)).
    
    Let a⋆ = s+t⋅α(s)/2 and b⋆=1-s+t⋅α(1-s)/2.
    We have formally that A = γ∑bᵗₙ/n^a⋆ = γ∑bᵗₙn^(b⋆-a⋆)/n^b⋆.
    Also formally, b⋆-a⋆ = 1-2s+t(α(1-s)-α(s))/2 = y-κ̄+ix+t⋅(α(s̄)-α(s))/2
    where by (18) pp 4, κ := t⋅(α(s)-α(1-s̄))/2 hence -κ̄ = t⋅(α(1-s)-α(s̄))/2.
    Thus we have
        bᵗₙ/n^a⋆ = {bᵗₙ⋅n^(y-κ̄)*n^(ix+(α(s̄)-α(s)))}/n^b⋆ = {bᵗₙ⋅n^(y-κ̄)*n^(i(x+ℑ(α(s̄)))}/n^b⋆
    where |n^(ix+ℑ(α(s̄))| = 1 since the exponent is purely imaginary.

    Note that, formally, αₙ = |γ|bᵗₙ⋅n^(y-κ̄) is the coefficient of the Dirichlet series in n^b⋆.
    """
function αterm(t::Real,n::Int, x::Real, y::Real)
    k = κ(t,x,y)
    s = s⁺(x,y)
    return bᵗₙ(t,n)*big(e)^(log(n)*(y-k')), big(e)^(log(n)*(im*+t*imag(α(s'))))
end

""" bound78(t::Real, x::Real, y::Real, λ::Vector{Real})

    """
function bound78(t::Real, x::Real, y::Real, λ̂::Vector{T}) where {T <: Real}
    N₀= N(t,x)
    s = s⁺(x,y)
    σ = real(s)
    α̂ = [αterm(t,n,x,y)[1] for n in 1:N₀]
    β̂ = [bᵗₙ(t,n) for n in 1:N₀]
    α̃ = Dirichlet_convolution(α̂,λ̂)
    β̃ = Dirichlet_convolution(β̂,λ̂)
    k = abs(κ(t,x,y))
    idx = find((α̃ .!= 0.0) .| (β̃ .!= 0.0))
    num = 1-real(α̃[1]) # ℑ(α̃[1]) is 0.0
    μ = num/(1+real(α̃[1]))
    for n in idx[2:end]
        num -= max(abs(β̃[n]-α̃[n]),μ*abs(β̃[n]+α̃[n]))/big(e)^(σ*log(n))
    end
    idx = find(λ̂ .!= 0)
    den = 0.0
    for d in idx
        den += abs(λ̂[d])/big(e)^(σ*log(d))
    end
    err = 0
    for n in 1:N₀
        err += bᵗₙ(t,n)*(big(e)^(k*log(n))-1)/big(e)^((σ-y)*log(n))
    end
    return num/den -abs(γₜ(t,x,y))*err 
end
                          
