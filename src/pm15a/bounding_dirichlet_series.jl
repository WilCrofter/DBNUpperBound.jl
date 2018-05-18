
export Dirichlet_convolution, ✪

function Dirichlet_convolution(f̂::Vector{T1}, ĝ::Vector{T2}) where {T1<:Number, T2<:Number}
    N=length(f̂)
    D=length(ĝ)
    ans = zeros(promote_type(T1,T2),N*D)
    for a in 1:N, b in 1:D
        ans[a*b] += f̂[a]*ĝ[b]
    end
    return ans
end

✪ = Dirichlet_convolution
