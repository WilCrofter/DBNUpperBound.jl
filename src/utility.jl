
function Φ_term{T<:AbstractFloat}(u::T,n::Int, PI)
    x = PI*n^2*exp(4*u)
    return (2*PI^2*n^4*exp(9*u-x)-3*PI*n^2*exp(5*u-x))
end

function Φ{T<:AbstractFloat}(u::T;n_max::Int=100, PI=T(π))
    running_sum = T(0)
    for n in 1:n_max
        running_sum += Φ_term(u,n,PI)
    end
    return running_sum
end

phi_decay = Φ #alias
KM_PI = 3.14159265358979323846
