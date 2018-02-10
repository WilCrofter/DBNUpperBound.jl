

function real_part_Γ(z)
    quadgk(t->exp(-t^2)*t^(2*real(z)-1)*cos(2*imag(z)*log(t)),0.0,Inf)
end

function imag_part_Γ(z)
    quadgk(t->exp(-t^2)*t^(2*real(z)-1)*sin(2*imag(z)*log(t)),0.0,Inf)
end

function test_Γ(z)
    ans=quadgk(t->exp(-t^2)*t^(2*z-1),0.0,Inf)
    2*ans[1],ans[2]
end
