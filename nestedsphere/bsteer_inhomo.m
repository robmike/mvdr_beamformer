function bsteer = bsteer_inhomo(soidir,numants)
pvar = 0.04; % Vorobyov (His variance is quoted in radians or degrees?)
theta = sqrt(pvar).*randn(numants,1);

bsteer = exp(j*theta).*bsteerdir(soidir,numants);

