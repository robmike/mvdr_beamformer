function ber = get_ber(sigarray, bitcombos,coeffs,noisevar)

% Analytically compute bit error rate

outvar=noisevar*coeffs'*coeffs; % Variance at output of filter

% countbits=numbits*numusers;
% 
 numbitcombo=size(bitcombos,1);

mmean = coeffs'*sigarray'*bitcombos; % Row vector (outputs of filter for each bit combo)
% ber = sum(Q(-real(mmean)/sqrt(outvar)))/numbitcombo;
ber = mean(Q(-real(mmean)/sqrt(outvar)));

