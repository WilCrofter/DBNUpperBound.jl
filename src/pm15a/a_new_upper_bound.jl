"""
    (82) pp 31
    """
function δ₁(t::Real, x::Real)
    in_region_5(t,x,.5) || error("t and x are not in region (5)")
    return ((t/4*log(x/(4*π)))^2 + 0.626)/(x-6.66)
end

"""

    (81) pp 31
    """
function F(N::Int, t::Real, σ::Real)
    # Σ bᵗₙ/n^σ = Σ exp(log(bᵗₙ)-σ*log(n))
    ans = 0.0
    for n in 1:N
        ans += big(e)^(logbᵗₙ(t,n)-σ*log(n))
    end
    return ans
end
