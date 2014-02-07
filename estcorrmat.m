function corrmat = estcorrmat(sigarray,numdatasamples,noisevar,varargin) 
% Estimate correlation matrix using sample mean with
% 'numdatasamples' data samples

delayspread=1;
numants=1;
numusers=size(sigarray,2);
spreadgain=size(sigarray,1);

if(nargin>3)
    sysparams = varargin{1};
    if(isfield(sysparams,'delayspread'))
        delayspread=sysparams.delayspread;
    end
    if isfield(sysparams,'numants')
        numants=sysparams.numants;
    end
    if isfield(sysparams,'numusers')
        numusers=sysparams.numusers;
    end
    if isfield(sysparams,'spreadgain')
        spreadgain=sysparams.spreadgain;
    end
end


if(delayspread==1)
    numbits=1; % Number of bits that affect the bit of interest
else
    numbits=3;
end
 

% - Generate array of random bits
randinfobits = 2*randbin(numusers*numbits,numdatasamples)-1;
% --------------------------------------

% Each column represents randomly generated noiseless received vectors
% and is of length numants*(delayspread+spreadgain-1)
x = sigarray*randinfobits;

numrows = numants*(delayspread+spreadgain-1);

x = x + sqrt(noisevar)*randn(numrows,numdatasamples);

% Sample mean of autocorrelation
corrmat = 1/numdatasamples * x*x';
