function [condcoeffs,ucoeffs,auxvecs] = ...
    auxvec_coeffs(corrmat,rakevec,numauxvecs);
% Returns the coefficients obtained using the auxiliary vector method and
% also the coefficents obtained with the unconditional auxiliary vector 
% method using the auxiliary vectors generated with the conditonal method
% the 


% Conditional opt
[condcoeffs,auxvecs] = condauxvec_coeffs(corrmat,rakevec,numauxvecs);

% Unconditional opt
ucoeffs = uncondauxvec_coeffs(corrmat,auxvecs,rakevec);