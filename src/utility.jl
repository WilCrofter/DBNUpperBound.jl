#= def phi_decay(u,n_max=100):
    running_sum=0
    for n in range(1,n_max+1):
        term1=2*PI_sq*pow(n,4)*exp(9*u) - 3*PI*pow(n,2)*exp(5*u)
        term2=exp(-1*PI*pow(n,2)*exp(4*u))
        running_sum += term1*term2
        #print n,term1, term2, running_sum
    return running_sum
=#

function Φ_term{T<:AbstractFloat}(u::T,n::Int)::T
    Tπ = T(π) # compute π to type T precision, e.g., BigFloat
    e5u = exp(5*u)
    return (2*(Tπ^2)*(n^4)*exp(9*u)-3*Tπ*(n^2)*e5u)*exp(-Tπ*n^2*e5u)
end

function Φ{T<:AbstractFloat}(u::T;n_max::Int=100)::T
    running_sum = T(0)
    for n in 1:n_max
        running_sum += Φ_term(u,n)
    end
    return running_sum
end

phi_decay = Φ #alias
