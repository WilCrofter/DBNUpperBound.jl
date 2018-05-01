#= Elementary estimates: random checks for Lemma 5.1 pp 11

The claims of Lemma 5.1 don't need computational checks, but a "computational version" of the proof of lemma_51vi does present cases in which provably correct O≤ relations fail computationally. 

Via [Elementary Asymptotics](http://michaelnielsen.org/polymath1/index.php?title=Effective_bounds_on_H_t_-_second_approach), "We use O≤(X) to denote any quantity bounded in magnitude by at most X. An equality A=B using this notation means that any quantity of the form A is also of the form B, though we do not require the converse to be true (thus we break the symmetry of the equality relation)." 

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
    # 1. Last step of the proof
    all(log.(x) .≥ a./b)  &&
        # 2. Even though the following is formally equivalent to 1,
        # it sometimes yields false, hence is one possible source of numerical trouble.
        #  all(b.*x.*log.(x) .≥ a.*x)
        # 3. Log equivalent of 1.
        all(log.(log.(x)) .≥ log.(a./b)) &&
        # 4. Adding log.(x) to both sides of 3 is formally equivalent to 2: xlog(x) ≥ xa/b 
        all(log.(x)+log.(log.(x)) .≥ log.(x) + log.(a./b)) &&
        # 5. Splitting out -log.(b) sometimes yields false
        # all(log.(x)+log.(log.(x)) .≥ log.(x) + log.(a)-log.(b))
        # 6. Replacing log(x) by log(x-c) on the RHS is formally equivalent to xlog(x) ≥ (x-c)a/b
        # hence to a/(xlog(x))- b/(x-c) ≤ 0.0 which is the essential step of the proof.
        all(log.(x)+log.(log.(x)) .≥ log.(x-c) + log.(a./b))
        # 7. the straightforward version of the formal inequality sometimes yields false.
        # all(a./(x.*log.(x)) .≤ b./(x-c))
end
    
