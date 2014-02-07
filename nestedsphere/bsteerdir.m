function beamsteer = bsteerdir(theta,numants)
% Generate beamsteering vectors steered to directions given by
% theta
% If theta is a vector return length(theta) beamsteering vectors
% Operation not defined in case theta a matrix
arrayspace = 1/2; % Array element spacing in wavelengths

[m,n] = size(theta);

%{
if(~isvector(theta))
    error('First argument must be a vector')
end
if(m<n) % Convert to column vector
    theta = theta';
end
%}

numusers = size(theta,1);
delayspread = size(theta,2);


% Angle of arrival of multipaths in radians. 
% Entry (i,j) is angle of jth multipath of ith user
%mpath_aoa = zeros(numusers,delayspread); % For no steering
% mpath_aoa=pi*rand(numusers,delayspread)-pi/2;
mpath_aoa = theta;

% Matrix multipath AOA augmented to 3 dimensions. third dimension is user,
% first dimension is antenna and third is multipath
mpath_aoa_aug = repmat(mpath_aoa,[1,1,numants]); % numsers-by-delayspread-by-numants
mpath_aoa_aug = permute(mpath_aoa_aug,[1 3 2]); % numusers-by-numants-by-delayspread
beamsteer = exp(  j*2*pi*arrayspace*...
    repmat(0:(numants-1),[numusers,1,delayspread]).* sin(mpath_aoa_aug)   ); 

% Change the indexing because we need to index by  mpath_aoa_aug(:,:,1) 
beamsteer = permute(beamsteer,[2 3 1]); % numants-by-delayspread-by-numusers

if(m==1 || n==1 || numants==1) % 
    beamsteer = squeeze(beamsteer);
end

% No steering can be achieved by:
%beamsteer=ones(numusers,numants,delayspread); % No steering
