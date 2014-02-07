function B = weightcols(weights,A)
% Weight the columns of the matrix A by the coefficients 'weights'

[m,n] = size(weights);

if m>n
   weights=weights';    
end

%{
if size(weights,2)~=1
    error('Weights must be a vector')
end
%}

[nrows,ncols]=size(A);
%B = repmat(weights,nrows,1).*A;
B = weights(ones(nrows,1),:).*A; % Faster than using repmat