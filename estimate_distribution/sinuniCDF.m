function z = sinuniCDF(y)
    
    [m,n] = size(y);
    z = zeros(m,n);
    z(-1<y & y<1) = asin(y(-1<y & y<1))/pi + 0.5;
    z(y>1) = 1;