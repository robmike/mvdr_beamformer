function sinr = get_sinr(coeffs,bsteer,soivar,intnocorrmat)
    sinr=...
        soivar*abs(coeffs'*bsteer)^2/real(coeffs'*intnocorrmat*coeffs);
end
