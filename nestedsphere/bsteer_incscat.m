function [outsteer] = bsteer_incscat(soidir,nants,npaths,nsamples)

std = 2;
m = soidir;
theta = std*randn(1,npaths-1) + m;
theta = theta*pi/180;

fvar = 1;
ovec = [1 npaths nsamples];
fad = sqrt(fvar).*complex(randn(ovec),randn(ovec));
fad = repmat(fad,[nants,1,1]);

steermat = repmat([bsteerdir(soidir,nants) bsteerdir(theta,nants)],[1,1,nsamples]);

outsteer = squeeze(sum( fad .*steermat,2));

%if nargout>1
% Return signal correlation matrix also
% sigcorrmat
%end