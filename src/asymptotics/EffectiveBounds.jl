
function eps_prime(t::Ty1, T::Ty1, s::Ty2, N::Int) where {Ty1<:Real, Ty2<:Number}
    α1 = αeff(s)
    real(s)+t*real(α1-t/4*log(N)) > 1 || error("argument s, ",s," is out of range.")
    ans = 0.0
    for n in 1:N
        α2 = (t^2/4)*abs(α1-log(n))^2 + 1/3 + t
        # For clarity:
        # ans += 1/n^(real(s) + t*real(α1)/2 - t/4*log(n)) * bigexp(α2/(2*(T-3.33))) * α2
        # For numerics:
        ans += bigexp(-n*(real(s) + t*real(α1)/2 - t/4*log(n)) + α2/(2*(T-3.33))) * α2
    end
    return real(ans)/2
end

function E1(t::Ty1, T::Ty1, s::Ty2, N::Int) where {Ty1<:Real, Ty2<:Number}
    return real(1/(8*(T-3.33))*bigexp(t/4*real(αeff(s))^2) * abs(H01(s))*eps_prime(t,T,s,N))
end

function E2(t::Ty1, T::Ty1, s::Ty2, N::Int) where {Ty1<:Real, Ty2<:Number}
    return E1(t,T,(1-s)',N)
end

function fσ(t::Ty1, σ::Ty2) where {Ty1<:Real, Ty2<:Real}
    y = 1-σ
    return 1/(2*√(π*t))*(bigexp(-(σ-(1-y)/2)^2/t) + bigexp(-(σ-(1+y)/2)^2)/t)
end

function vσ(t::Ty1, T0::Ty2, σ::Ty3) where {Ty1 <: Real, Ty2 <: Real, Ty3 <:Real}
    T0prime = T0 + π*t/8
    a0 = √(T0prime/(2*π))
    if σ ≥ 0
        return 1 + (0.400*9^σ/a0 + 0.346*2^(3σ/2)/a0^2)
    else
        psum = 0.0
        for k in 1:1:floor(Int,(4-σ))
            psum += (1.1)^k*bigΓ(k/2)/a0^k
        end
        return (9/10)^ceil(Int,-σ)*psum
    end
end

function wσ(t::Ty1, T0::Ty2, σ::Ty3) where {Ty1 <: Real, Ty2 <: Real, Ty3 <:Real}
    T0prime = T0 + π*t/8
    a0 = √(T0prime/(2*π))
    xpo = max(0,(σ-1))/4*log(1 + σ^2/T0prime^2) + (σ<0 ? T0prime/2*atan(σ/T0prime) : 0) + 1/(12*(T0prime-0.33))
    return real((1+σ^2/T0prime^2)^(1/2)*(1 + (1-σ)^2/(T0prime^2))^(1/2)*bigexp(xpo))
end

function E3(t::Ty1, T::Ty1, T0::Ty1, s::Ty2; limit=100.0) where {Ty1 <: Real, Ty2 <: Number}
    T0 ≥ 10 || error("T0 = ",T0," is out of applicable range.")
    T ≥ T0 || error("T = ", T, " must be at least ", T0, " (T0).")
    Tprime = imag(s) + π*t/8
    t ≤ 1/2 || error("t must not exceed 1/2. Given value is ",t,".")
    σ = real(s)
    # quadgk will balk at limits of integration ≥ 100.
    defint = quadgk((σ)->vσ(t,T0,σ)*wσ(t,T0,σ)*fσ(t,σ),-limit,limit)
    μ = real(1/8*√(π)bigexp(-t*π^2/64-π*T/4)*Tprime^3/2)
    return μ*real(defint[1]),μ*defint[2]
end

