module Numerics

function myΓ(z)
    return quadgk((x)->x^(z-1)*exp(-x), 0.0, Inf, abstol=eps(abs(z)))
end

end
