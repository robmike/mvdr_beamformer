function [w,variance]=mvdr_coeff(R,s)
% Returns the coefficients for the MVDR (minimum variance distortionless 
% response filter) for the signature vector s and correlation matrix R

Rinv=inv(R);
w=(Rinv*s)/(s'*Rinv*s);

if nargout>1
   variance=(s'*Rinv*s)^(-1); % Variance at output of MVDR filter
end

