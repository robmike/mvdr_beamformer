function [coeffs,auxvecs]=condauxvec_coeffs(R,steer,numvecs)
% AUXVEC_COEFFS_COND
% Returns filter coefficients for the conditionally optimized auxiliary 
% vector method as well as the computed auxiliary vectors
% numvecs - Number of auxiliary vectors

steer = steer/norm(steer); % Rake steering vector needs to be normalized

numrows=size(R,2);

weights=zeros(1,numvecs);
auxvecs=zeros(numrows,numvecs);

for i=1:numvecs
[weights(i) auxvecs(:,i)] = nextparams(weights,R,auxvecs,steer);
end

coeffs = steer - sum( weightcols(weights,auxvecs) , 2);


end


function [ci,auxveci]=nextparams(weights,R,auxvecs,steer)
% Return the ith filter coefficent and ith auxiliary vector 
% given coefficients 1,...,i-1 and auxiliary vectors 1,...,i-1

auxveci = nextauxvec(weights,R,auxvecs,steer); % ith auxiliary vector

wsum = sum( weightcols(weights,auxvecs) , 2);
wsum = steer - wsum;
ci = auxveci' * R * wsum;
ci = ci/(auxveci'*R*auxveci);


end