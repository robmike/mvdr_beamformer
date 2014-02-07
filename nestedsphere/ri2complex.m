function y = ri2complex(x)
% Convert real vector of dimension 2*n which represents a n-dimensionsal
% complex into said vector.
% First n entries of x are real parts, last n are complex parts
y = complex(x(1:end/2,:),x(end/2+1:end,:));
end