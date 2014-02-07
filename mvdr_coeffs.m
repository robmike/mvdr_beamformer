function [coeffs,output_variance]=mvdr_coeffs(R,s)
% MVDR_COEFFS
% Returns the coefficients for the MVDR (minimum variance distortionless 
% response filter) for the steering vector s and correlation matrix R

% Rinv=inv(R);
% coeffs=(Rinv*s)/(s'*Rinv*s);
% 
% if nargout>1
%    output_variance=(s'*Rinv*s)^(-1); % Variance at output of MVDR filter
% end

w=R\s;
coeffs=(w)/(s'*w);

if nargout>1
   output_variance=(s'*w)^(-1); % Variance at output of MVDR filter
end

