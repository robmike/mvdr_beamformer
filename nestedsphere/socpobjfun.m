function outsinr = socpobjfun(x,bsteerest_ri,R,tbsteer,soivar,intnocorrmat)
% Returns output SINR for a given

radius = torow(x(1:end/2));
atten = torow(x(end/2+1:end));
y1 = socp_coeffs(bsteerest_ri,atten,R,radius);
y1 = complex(y1(1:end/2),y1(end/2+1:end));
outsinr = get_sinr(y1,tbsteer,soivar,intnocorrmat);
