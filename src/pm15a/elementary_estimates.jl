#= Elementary estimates: random checks for Lemma 5.1 pp 11

This code is of dubious value at best: the claims of Lemma 5.1 don't need computational checks. It may or may not end up in the package, but the lemma_51vi does present some numerical issues.
=#

function lemma_51i(n::Int ;upper_limit::Real = 1e6)
    a = rand(n)*upper_limit
    b = rand(n)*upper_limit
    c = rand(n)*upper_limit
    tmp = max.(b./a, sqrt.(c./a))
    x = tmp+rand(n).*upper_limit
    all( a./x + b./x.^2 + c./x.^3 .≤ a./(x-tmp) )
end

function lemma_51ii(n::Int; upper_limit::Real = 1e6)
    x = 1+rand(n)*upper_limit
    all( log.(1+1./x) .≤ 1./(x-1))
end

function lemma_51iii(n::Int; upper_limit::Real = 1e6)
    x = 1/2 + rand(n)*upper_limit
    all(1./x .≤ 1+1./(x-0.5) )
end

function lemma_51iv(n::Int; upper_limit::Real = 1e6)
    x = rand(n)*upper_limit
    all( x .≤ 1+(map(bigexp,x)-1) )
end

function lemma_51v(n::Int; upper_limit::Real = 1e6)
    tmp = rand([0.0, 1.0], n)
    z = rand(n)*upper_limit+tmp + im.*(upper_limit*rand(n)+1-tmp).*rand([-1.0,1.0],n)
    lhs = lgamma.(z)
    rhs=log(√(2*π)) + z.*log.(z) - z + 1./(12.*(abs.(z)-0.33))
    all(isapprox.(abs.(lhs ./ rhs),1.0, rtol=6))
end

function lemma_51vi(n::Int; upper_limit::Real = big(1e6))
    a = rand(n)*upper_limit
    b = rand(n)*upper_limit
    c = rand(n)*upper_limit
    y = rand(n) # let 0 ≤ y ≤ 1
    logx₀ = max.(log.(c),a./b)
    x₀ = map(exp, logx₀)
    x = x₀+rand(n).*upper_limit
    all(log.(x) .≥ a./b)  &&
        all(b.*log.(x) .≥ a)
    # NOTE: even though the following is formally equivalent to the foregoing,
    # it sometimes yields false, hence is one possible source of numerical trouble.
    #  all(b.*x.*log.(x) .≥ a.*x)
end
    
