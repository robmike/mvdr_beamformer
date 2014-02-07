function [outsinr radius atten coeffs] = optsocp(numcon,bsteerest_ri,R,tbsteer,soivar,intnocorrmat)
% Returns the optimum socp coefficients for a given SNR and data
% correlation.
n = numcon;
M = length(bsteerest_ri)/2;
objfunc = @(x) -socpobjfun(x,bsteerest_ri,R,tbsteer,soivar,intnocorrmat);
x0 = [1:numcon,100.^(0:-1:-numcon+1)]';
reps = 1e-3; % Separation
A = zeros(2*n-2,2*n);
C = [diag(ones(n-1,1)),zeros(n-1,1)]+[zeros(n-1,1),diag(-ones(n-1,1))];
A(1:n-1,1:n) = C;
A(n:end,n+1:end) = -C;
b = [-reps(ones(n-1,1)); -reps(ones(n-1,1))];
lbound = [repmat(0.1,n,1); repmat(10^-8,n,1)];
ubound = [sqrt(M(ones(n,1))) - 0.5; ones(n,1)];
options = optimset('display','off');
[xhat outsinr] = fmincon(objfunc,x0,A,b,[],[],lbound,ubound,[],options);
outsinr = -outsinr; % Minimized negative of objective function
radius = torow(xhat(1:end/2));
atten = torow(xhat(end/2+1:end));

if nargout>3
   coeffs = socp_coeffs(bsteerest_ri,atten,R,radius); 
end