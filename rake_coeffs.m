function coeffs=rake_coeffs(channel,energy,signature,beamsteer)
% channel - Vector of channel coefficients for user of interest, one entry per multipath 
% energy - Transmitted energy of user of interest
% signature - signature waveform of user of interest, sampled at chip rate. Must have unit energy. Must account for any necessary time delays (e.g. asynchronisms)
% beamsteer - Array of beamsteering vectors as columns. nth col is beamsteering vector for nth multipath
% Signature waveform is assumed to be bandlimited to 1/chip_period

% Todo: make everything the expected dimension (row or col vector)

% npaths=size(channel,1); % Number of multipaths for user of interest

% For each antenna find response of signature to channel
% i.e. convolve each row of beamsteer*channel with signature
coeffs = sigout(channel,energy,signature,beamsteer);
coeffs = coeffs(:); % Stack each column on top of each other and into one big vector


% Comment this out because we don't want to specify correlation matrix R as input (it is unnecessary for calculating coefficients)
%if nargout>1
%    output_var = coeffs' * R * coeffs; % Variance at output of rake filter
%end

