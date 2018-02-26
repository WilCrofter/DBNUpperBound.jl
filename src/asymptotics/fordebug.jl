""" Ht_AFE_B(z, t)
    For comparison, line-by-line port of python function of the same name.
"""

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
    return B_pre*B_sum
end

