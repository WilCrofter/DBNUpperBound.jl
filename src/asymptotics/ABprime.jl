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


function Ht_AFE_B(z, t)
    littles = (1+im*real(z)-imag(z))/2
    z = big(real(z))+im*(big(imag(z)))
    t = big(t)
    s = (1+im*real(z)-imag(z))/2
    tau = sqrt(imag(s)/(2*big(π)))
    M = floor(Int,tau)
    B_pre = (1/16) * s * (s-1) * big(π)^((s-1)/2) * bigΓ((1-littles)/2)
    B_sum = 0.0
    for m in 1:M
        B_sum += bigexp((t/16) * log((5-s)/(2 * big(π) * m^2))^2)/m^(1-s)
    end
    B_pre*B_sum
end
#=
x
	B0
	B′0
	Beff0
	Btoy0
103
	(3.4405+3.5443im)*1e−167
	(3.4204+3.5383im)*1e−167
	(3.4426+3.5411im)*1e−167
	(2.3040+2.3606im)*1e−167
104
	(−1.1843−7.7882im)*1e−1700
	(−1.1180−7.7888im)*1e−1700
	(−1.1185−7.7879im)*1e−1700
	(−1.1155−7.5753im)*1e−1700
105
	(−7.6133+2.5065im)∗1e−17047
	(−7.6134+2.5060im)∗1e−17047
	(−7.6134+2.5059im)∗1e−17047
	(−7.5483+2.4848im)∗1e−17047
106
	(−3.1615−7.7093im)∗1e−170537
	(−3.1676−7.7063im)∗1e−170537
	(−3.1646−7.7079im)∗1e−170537
	(−3.1590−7.6898im)∗1e−170537
107
	(2.1676−9.6330im)∗1e−1705458
	(2.1711−9.6236im)∗1e−1705458
	(2.1571−9.6329im)∗1e−1705458
	(2.2566−9.6000im)∗1e−1705458
=#
