function beamsteer = bsteeruni(numusers,delayspread,numants)
% Uniform AOA for beamsteering vector

arrayspace = 1/2; % Array element spacing in wavelengths
% Angle of arrival of multipaths in radians. 
% Entry (i,j) is angle of jth multipath of ith user
%mpath_aoa = zeros(numusers,delayspread); % For no steering
mpath_aoa=pi*rand(numusers,delayspread)-pi/2;

% Matrix multipath AOA augmented to 3 dimensions. third dimension is user,
% first dimension is antenna and third is multipath
mpath_aoa_aug = repmat(mpath_aoa,[1,1,numants]); % numsers-by-delayspread-by-numants
mpath_aoa_aug = permute(mpath_aoa_aug,[1 3 2]); % numusers-by-numants-by-delayspread
beamsteer = exp(  j*2*pi*arrayspace*...
    repmat(0:(numants-1),[numusers,1,delayspread]).* sin(mpath_aoa_aug)   ); 

% Change the indexing because we need to index by  mpath_aoa_aug(:,:,1) 
beamsteer = permute(beamsteer,[2 3 1]); % numants-by-delayspread-by-numusers


% No steering can be achieved by:
%beamsteer=ones(numusers,numants,delayspread); % No steering