#= Note that we generally assume s = σ+im*T = (1+im*(x+im*y))/2. Since it is not entirely clear at the the time of this writing, that this is always the case some alternative forms are provided for ambiguous cases.
=#

function eA(t::Real, x::Real, y::Real)
    s = (1+im*(x+im*y))/2
    ans = 0.0
    for n in 1:N(t,x)
        ans += n^y*bᵗₙ(t,n)/n^(real(s)+real(κ(t,x,y)))*ϵₜₙ((1-y+im*x)/2)
    end
    return abs(γ(t,x,y))*ans
end

function eA(t::Real, z::Number)
    return eA(t, real(z), imag(z))
end

function eB(t::Real, x::Real, y::Real)
    s = (1+im*(x+im*y))/2
    ans = 0.0
    for n in 1:N(t,x)
        ans += bᵗₙ(t,n)/n^real(s)*ϵᵗₙ((1+y+im*x)/2)
    end
    return ans
end

function eB(t::Real, z::Number)
    return eB(t,real(z),imag(z))
end

function eC(t::Real, x::Real, y::Real, σ::Real, T::Real)
    T′ = T - π*t/8
    return exp(t*π^2/64)*abs(M₀(im*T'))/abs(Mₜ((1+y+im*x)/2))*(ϵ̃((1-y+im*x)/2) + ϵ̃((1+y+im*x)/2))
end

function eC(t::Real, z::Number, σ::Real, T::Real)
    return eC(t,real(z),imag(z),T)
end

function eC(t::Real, x::Real, y::Real, s::Number)
    return eC(t,x,y,real(s),imag(s))
end

function eC(t::Real, z::Number, s::Number)
    return eC(t,real(z),imag(z),real(s),imag(s))
end

function eC(t::Real, x::Real, y::Real)
    s = (1+im*(x+im*y))/2
    return eC(t,x,y,real(s),imag(s))
end

function eC(t::Real, z::Number)
    s = (1+im*z)/2
    return eC(t,real(z),imag(z),real(s),imag(s))
end

function eC0(t::Real, x::Real, y::Real, T::Real)
    T′ = T - π*t/8
    return exp(t*π^2/64)*abs(M₀(im*T'))/abs(Mₜ((1+y+im*x)/2))*(1 + ϵ̃((1-y+im*x)/2) + ϵ̃((1+y+im*x)/2))
end

function eC0(t::Real, z::Number, T::Real)
    return eC0(t,real(z),imag(z),T)
end

function eC0(t::Real, x::Real, y::Real)
    T = imag((1+im*(x+im*y))/2)
    return eC0(t,x,y,T)
end

function eC0(t::Real, z::Number)
    T = imag((1+im*z)/2)
    return eC0(t,real(z),imag(z),T)
end
