export δ₁, δ₁_ub, bound23, bound23a, util_23b, bound23b, bound85

""" δ₁_ub(t₀::Real, x₀::Real)

    Returns upper bound (83) pp 31 for δ₁ valid for 0≤t≤t₀ and x≥x₀ in region 5. Since, in region 5, δ₁ is monotone decreasing in x, and monotone increasing in t, the value of δ₁ at t₀ and x₀ serves.

    (83) pp 31
    """
function δ₁_ub(t₀::Real, x₀::Real)
    return δ₁(t₀,x₀)
end


""" δ₁(t::Real, x::Real)

    Returns δ₁ := ((t/4⋅log(x/(4π)))² + 0.626)/(x-6.66) as defined in (82) pp 31
    """
function δ₁(t::Real, x::Real)
    (0≤t≤0.5)&&(x≥200.0) || error("t and x are not in region (5)")
    return ((t/4*log(x/(4*π)))^2 + 0.626)/(x-6.66)
end

""" bound23(t::Real, x::Real, y::Real)

    Returns an upper bound to eA+eB as per inequality (23) pp. 4
        eA+eB ≤ ∑ (1+|γ|N^|κ|n^y)bᵗₙ/n^(ℜ(s⋆)))⋅(e^δ₁-1)
    """
function bound23(t::Real, x::Real, y::Real)
    s = s⁺(x,y)
    star = real(1-s + t/2*α(1-s))
    absγ = abs(γₜ(t,x,y))
    absκ = abs(κ(t,x,y))
    N₀ = N(t,x)
    μ = N₀^(absκ)*absγ 
    ans = 0.0
    for n in 1:N₀
        ans += (1+μ*big(e)^(y*log(n)))*big(e)^(logbᵗₙ(t,n)-star*log(n))*(big(e)^δ₁(t,x)-1)
    end
    return ans
end

""" bound23a(t₀::Real, x₀::Real, t::Real, x::Real, y::Real)

    A cruder form of bound23 based on nʸ≤Nʸ when y≥0, 1≤n≤N, and on an upper bound for δ₁ given by (83) pp. 31:
   
    ∑ (1+|γ|N^|κ|n^y)bᵗₙ/n^(ℜ(s⋆)))⋅(e^δ₁-1) ≤ (1+|γ||N|^(|κ|+y))⋅(e^δ₁(t₀,x₀) -1)⋅∑ bᵗₙ/n^(ℜ(s⋆))
    """
function bound23a(t₀::Real, x₀::Real, t::Real, x::Real, y::Real)
    (0≤t≤t₀) && (x≥x₀) && (y≥0) || error("t or x is out of range; 0≤t≤t₀, x≥x₀, and y≥0 are required.")
    s = s⁺(x,y)
    absγ = abs(γₜ(t,x,y))
    absκ = abs(κ(t,x,y))
    N₀ = N(t,x)
    return (1+absγ*big(e)^((absκ+y)*log(N₀)))*(big(e)^δ₁_ub(t₀,x₀)-1)*F(N₀,t,real(1-s + t/2*α(1-s)))
end


"""  F(N::Int, t::Real, σ::Real)

    Returns Σ bᵗₙ/n^σ = Σ exp(log(bᵗₙ)-σ*log(n) as per (81) pp 31
    """
function F(N::Int, t::Real, σ::Real)
    ans = 0.0
    for n in 1:N
        ans += big(e)^(logbᵗₙ(t,n)-σ*log(n))
    end
    return ans
end

""" util_23b(t::Real, x::Real, y::Real)

    Returns exp(0.02*y + ty⋅(log(x/(4π))+ log(1+πt/(4x))/(4(x-6)) which is an upper bound to |γ||N|^(|κ|+y)⋅(e^δ₁(t₀,x₀)-1) and which is decreasing  x (x≥e), increasing in y and t, hence valid for x′≥x, 0≤y′≤y, and 0≤t′≤t.

    It is assumed that N := √(x/(4*π) + t/16).
    """
function util_23b(t::Real, x::Real, y::Real)
    return big(e)^(0.02*y + t*y*(log(x/(4*π))+log(1+π*t/(4*x)))/(4*(x-6)))
end

""" bound23b(t₀::Real, x₀::Real, t::Real, x::Real, y::Real; y₀::Real=1)
    
    Returns (1+exp(0.02y₀ + t₀y₀⋅(log(x₀/(4π))+ log(1+πt₀/(4x₀))/(4(x₀-6)))))⋅Σ bᵗₙ/n^ℜ(s⋆), an upper bound to 23a and hence an upper bound to eA+eB.
    """
function bound23b(t₀::Real, x₀::Real, t::Real, x::Real, y::Real; y₀::Real=1)
    s=s⁺(x,y)
    return (1+util_23b(t₀,x₀,y₀))*F(N(t,x),t,real(1-s+t/2*α(1-s)))
end

""" bound85(t::Real, x₀::Real, y₀::Real, N₀::Int)

    Returns an upper bound to eC₀ as per inequality (85) pp. 32, valid for x≥x₀≥200, y≥y₀>0, and N≥N₀. Note that the bound is decreasing in t, hence valid for t′≥t.
    """
function bound85(t::Real, x₀::Real, y₀::Real, N₀::Int)
    (t>0)&&(x₀≥200) && (y₀>0) && (N₀≥1) || error("Parameters violate at least one of the conditions t>0, x₀≥200, y₀>0, N₀≥1.")
    l4 = log(x₀/(4*π))
    return big(e)^(-(1+y₀)/4*l4 - t/16*l4^2 +(3*abs(l4+im*π/2)+3.58)/(x₀-8.52)) * 
        (1+1.24*(3^y₀+3^(-y₀))/(N₀-1.25)+6.72/(x₀-6.66))
end
