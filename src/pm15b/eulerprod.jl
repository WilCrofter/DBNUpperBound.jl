
export eulerprod, adhoc
using Primes

""" eulerprod(x::Real, primes::Vector{T}=[i for i in 2:27 if isprime(i)]) where {T<:Integer}

    Return |Π(1/(1-p^(ix/2-1)))| where the product is over p∈primes (default primes ≤ 27). See pp. 32.
    """
function eulerprod(x::Real, primes::Vector{T}=[i for i in 2:n if isprime(i)]) where {T<:Integer}
    tmp = im*x/2-1
    return abs(prod(1 ./ (1-primes.^tmp)))
end

""" adhoc(;X::T=6e10, n::I=27, Q::I=10^5, set::Vector{T}=[-0.5,0.0,0.5], threshold::Real=4.0) where {T<:Real, I<:Integer}

    Execute the ad hoc procedure described on pp. 32-33 of the reference to determine placement of the barrier, returning accepted values of q and corresponding values of |f₀| as a 2-tuple. 

    X (default 6x10¹⁰) approximate location of the barrier.
    n (default 27) upper limit of primes to include in the euler product
    Q (defualt 10⁵) upper limit of natural numbers, q, in the expression x-X-q
    set (default [-0.5,0.0,0.5]) array of allowable values for x-X-q q∈{1,…,Q}
    threshold (default 3.944) acceptance threshold for values of eulerproduct(x-X-q)

    TIMING (2 core, 3.6 GHz, evaluation of f₀ accounting for most.) 
    julia> @time q,f=M.adhoc()
    275.663969 seconds (846.60 M allocations: 27.724 GiB, 10.44% gc time)

    NOTE: The reference reports that a threshold of 4.0 returns the 7 integers 1046, 22402, 24198, 52806,
77752, 83952, and 99108. With a threshold of 4.0 the current function returns only the first 5 of these. With the default threshold, 3.944, this function returns 11 which include the reported 7. The heuristically optimal value of q, that which maximizes the minimum value of |f₀(X+q+r+i)| where r ∈ set and i=√(-1), is the same however.

     q       f₀
  1046  4.24874
 22402  4.2109 
 24198  4.28588
 52806  4.11698
 57701  3.68531
 61762  4.18553
 72114  3.9124 
 72857  3.444  
 77752  3.93115
 83952  4.31925 *
 99108  3.80751

    """
function adhoc(;X::T=6e10, n::I=27, Q::I=10^5, set::Vector{T}=[-0.5,0.0,0.5], threshold::Real=3.944) where {T<:Real, I<:Integer}
    accepted = spzeros(Bool,Q)
    tmp = Vector{T}(length(set))
    primes =[i for i in 2:n if isprime(i)]   
    for q in 1:Q
        for i in eachindex(set)
            tmp[i]= eulerprod(X+q+set[i],primes)
        end
        if minimum(tmp) ≥ threshold
            accepted[q]=true
        end
    end
    accepted=find(accepted)
    f₀=Vector{typeof(X)}(length(accepted))
    for i in eachindex(accepted)
        for j in eachindex(set)
            tmp[j] = abs(fₜ(0,X+accepted[i]+set[j]+im)[1])
        end
        f₀[i]=minimum(tmp)
    end
    return accepted, f₀
end
    
