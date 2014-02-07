function vec=nextauxvec(weights,R,auxvecs,steer)
% Returns the 'next' auxiliary vector for the conditionally optimized
% auxiliary vector method
% Requires past auxiliary vectors and filter coefficients as well as
% steering vector


% Make weights a row vector
[nrows,ncols]=size(weights);
if nrows>ncols
    weights=weights';
end

[nrows,ncols]=size(auxvecs);

wsum = sum( weightcols(weights,auxvecs) , 2);
wsum = steer - wsum;

vec = R*wsum;
temp = steer'*R*wsum*steer;
vec = vec - temp;
temp = sum(weightcols(auxvecs'*(R*wsum),auxvecs),2);
vec = vec - temp;

vec=vec/norm(vec); % Normalize

end