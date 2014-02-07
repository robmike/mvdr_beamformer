function out = sigout(channel,energy,signature,beamsteer)

% Channel response of given signature waveform, sampled at chip rate
% Output is a matrix with (i,j) entry being the jth sample of the output at
% the ith antenna

%beamsteer=squeeze(beamsteer); % If beamsteer is 3-D make it 2-D
%out = conv2(signature,weightcols(sqrt(energy)*channel,beamsteer));

% Inline weightcols 
weights=sqrt(energy)*channel; % to weight each column of beamsteer
numants=size(beamsteer,1);
out = conv2(signature,weights(ones(numants,1),:).*beamsteer);