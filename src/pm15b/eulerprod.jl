
export eulerprod, adhoc
using Primes

""" eulerprod(x::Real, primes::Vector{T}=[i for i in 2:27 if isprime(i)]) where {T<:Integer}

    Return |Π(1/(1-p^(ix/2-1)))| where the product is over p∈primes (default primes ≤ 27). See pp. 32.
    """
function eulerprod(x::Real, primes::Vector{T}=[i for i in 2:n if isprime(i)]) where {T<:Integer}
    tmp = im*x/2-1
    return abs(prod(1./(1-primes.^tmp)))
end

""" adhoc(;X::T=6e10, n::I=27, Q::I=10^5, set::Vector{T}=[-0.5,0.0,0.5], threshold::Real=4.0) where {T<:Real, I<:Integer}

    Execute the ad hoc procedure described on pp. 32-33 of the reference to determine placement of the barrier.

    X (default 6x10¹⁰) approximate location of the barrier.
    n (default 27) upper limit of primes to include in the euler product
    Q (defualt 10⁵) upper limit of natural numbers, q, in the expression x-X-q
    set (default [-0.5,0.0,0.5]) array of allowable values for x-X-q q∈{1,…,Q}
    threshold (default 3.944) acceptance threshold for values of eulerproduct(x-X-q)

    NOTE: The reference reports that a threshold of 4.0 returns the 7 integers 1046, 22402, 24198, 52806,
77752, 83952, and 99108. With a threshold of 4.0 this function returns only the first 5 of these. With the default threshold, this function returns 11 which include the reported 7.
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
    return find(accepted)
end
    
