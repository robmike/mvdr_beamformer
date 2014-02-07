function z = normfunc2(x,theta,numants)
    
    M = numants;
    ed = bcast(x,-pi*sin(theta));
    z = sqrt(2*M-1-(2*M-1)*diric(ed,(2*M-1)));