function channel = raylchan(m,n,freq,phase,chanvar)
% m by n Rayleigh channel (zero mean complex Gaussian)
% with mx1 phase offset

channel = sqrt(chanvar)*...
    complex(randn(m,n),randn(m,n)).*...
    exp(-j* (2*pi*freq*repmat(1:n,m,1) + ...
    repmat(phase,1,n) ) );


% raylchan(numusers,delayspread,freq,totcarphase,chanvar)