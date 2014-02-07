function corrmat=idealcorr(sigarray,noisevar,varargin)
% Retrieve theoretical correlation matrix of received data samples

numusers=size(sigarray,2);

userrange = 1:numusers;

% This should be handled okay if varargin{1} is an empty matrix
if ~isempty(varargin)
    userrange=varargin{1}; % Userful for computing Interference correlation
    numusers=length(userrange);
end

S = sigarray(:,userrange);

% Bits from different users are uncorrelated as are bits in
% different bit intervals so we consider each bit interval for each
% user separately and add them.

corrmat =  S*S';


% Noise correlation matrix (white noise)
corrmat = corrmat + noisevar*eye(size(corrmat,2));
