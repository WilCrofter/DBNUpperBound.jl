# http://michaelnielsen.org/polymath1/index.php?title=Polymath15_test_problem

function utilA1(s::T) where {T <: Number}
    return (2/8)*sqrt(2*π)*bigexp(-(s/2)*log(π) + ((s+4/2)-1/2)*log((s+4)/2)-(s+4)/2)
end

function utilB1(s::T) where {T <: Number}
    return utilA1(1-s)
end

function utilA2(s::T1, t::T2, n::Int) where {T1 <: Number, T2 <: Real}
    return bigexp(-s*log(n) + (t/16)*log((s+4)/(2*π*n^2))^2)
end

function utilB2(s::T1, t::T2, n::Int) where {T1 <:Number, T2 <: Real}
    return utilA2(1-s, t, n)
end

function Aprime(t::T1,N::I1,s::NBC) where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    psum = 0
    for n in N:-1:1
        psum += utilA2(s,t,n)
    end
    return utilA1(s)*psum
end

function Bprime(t::T1,M::I1,s::NBC) where {T1<:Number,I1<:Integer,NBC<:NotBigComplex}
    psum = 0
    for m in M:-1:1
        psum += utilB2(s,t,m)
    end
    return utilB1(s)*psum
end
