function [cineq,ceq]=socpcon(z,c,At,K)
% function [c,ceq]=socpcon(z,c,At,K)
% Second order cone constraint function
% Too complicated to explain how it works -- No documentation for
% now :(


if(any(K.q<2) || any(fix(K.q)~=K.q))
    error('SOCPCON: Cone dimensions must be integers larger than 1');
end

if(~isfield(K,'f')) % Rows of c- At*x restricted to zero
    K.f = 0;
    ceq = 0;
elseif(~isscalar(K.f) || K.f<0 ||fix(K.f)~=K.f) 
    error('SOCPCON: K.f must be an integer greater than zero');
else
    ceq = zeros(K.f,1);
    
    % Equality constraints
    ceq(1:K.f) = c(1:K.f) - At(1:K.f,:)*z;
end

cineq = zeros(length(z),1);

sidx = 1;
for i = 1:length(K.q)
    n = K.q(i);
    coneh = c(sidx)-At(sidx,:)*z; % Cone height term
    rowrange = sidx+1:sidx+n-1;
    conenorm = norm(c(rowrange) - At(rowrange,:)*z); % Cone norm term
    cineq(i)= conenorm - coneh;
    sidx = sidx + n; % First row for this cone
end


