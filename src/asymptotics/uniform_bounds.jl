αeff = DBNUpperBound.Asymptotics.αeff # not exported

""" E2_by_B0(N,t,y)
    
    Returns a uniform upper bound for E2/B0eff on the interval xN(N) ≤ ℜ(z) < xN(N+1).

     For more documentation see notebooks/uniform_bounds.ipynb and the derivation, http://michaelnielsen.org/polymath1/index.php?title=Controlling_H_t-A-B/B_0#Estimation_of_E_1.2CE_2.
    """
function E2_by_B0(N::Int, t::Real, y::Real)
    y ≥ 1/3 || error("Bound applies only if y ≥ 1/3; y = ",y," is out of range.")
    x(n) = 4*π*n^2 - π*t/4 # note: a function
    s⁺N = (1+y+im*(x(N)+x(N+1))/2)/2
    κ = (x(N+1)-x(N))/(4*(x(N)-6))
    a(n) = abs(αeff(s⁺N)-log(n))
    c⁺(n) = t^2/4*(a(n)+κ)^2+1/3+t
    b = (1+y)/2 + t/2*real(αeff(s⁺N))-t/2*κ # a constant for the following sum
    psum = 0.0
    # Note the following expressions in T are maximized when T=x/2 is minimized hence
    T = x(N)/2
    for n in 1:N
        c = c⁺(n)
        psum += c*exp(c/(2*(T-3.33)))/n^(b-t/4*log(n)) 
    end
    return 1/(T-3.33)*psum/2
end

""" E1_by_B0(N,t,y)
    
    Returns a uniform upper bound for E1/B0eff on the interval xN(N) ≤ ℜ(z) < xN(N+1).

     For more documentation see notebooks/uniform_bounds.ipynb and the derivation, http://michaelnielsen.org/polymath1/index.php?title=Controlling_H_t-A-B/B_0#Estimation_of_E_1.2CE_2.
    """
function E1_by_B0(N::Int, t::Real, y::Real)
    y ≥ 1/3 || error("Bound applies only if y ≥ 1/3; y = ",y," is out of range.")
    x(n) = 4*π*n^2 - π*t/4 # note: a function
    s⁻N = (1-y+im*(x(N)+x(N+1))/2)/2
    κ = (x(N+1)-x(N))/(4*(x(N)-6))
    a(n) = abs(αeff(s⁻N)-log(n))
    c⁻(n) = t^2/4*(a(n)+κ)^2+1/3+t
    δ = π*y/(2(x(N)-6-(14+2*y)/π)) + 2*y*(7+y)/x(N)^2*log(abs(1+y+im*x(N+1))/(4*π))
    b = (1-y)/2 + t/2*real(αeff(s⁻N))-t/2*κ # a constant for the following sum
    psum = 0.0
    # Note the following expressions in T are maximized when T=x/2 is minimized hence
    T = x(N)/2
    for n in 1:N
        c = c⁻(n)
        psum += c*exp(c/(2*(T-3.33)))/n^(b-t/4*log(n)) 
    end
    return 1/(T-3.33)*exp(δ)*N^(-y)*psum/2
end


""" E3star_by_B0(N,t,y)

    Returns a 2-tuple representing E3*/E3main and E3main/|B₀eff| where E3main = (T0′)^(3/2)*exp(-πT0/4). Both quantities have been shown to be monotonically decreasing for N ≥ 3.

    For more documentation see notebooks/uniform_bounds.ipynb and the derivation, http://michaelnielsen.org/polymath1/index.php?title=Controlling_H_t-A-B/B_0#Estimation_of_E_1.2CE_2.
    """
function E3star_by_B0(N::Int, t::Real, y::Real)
    xN = 4*π*N^2 - π*t/4
    T0 = xN/2
    s = (1-y)/2 + im*T0
    T0′ = T0 + π*t/8
    a0 = √(T0′/(2*π))
    E3star_over_main = 1/8*√(π)*exp(-t*π^2/64 + 0.181/(T0′-3.33))*(1+5.15/(a0-1.25))
    main_over_B0 = (T0′)^(3/2)*real(bigexp(-π*T0/4))/abs(Beff(t,s,1))
    return E3star_over_main, main_over_B0
end

