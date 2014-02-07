function outsteer = bsteer_coscat(losdir,nants,npaths,nsamples)
% Simulate scattered beamsteer

npaths = npaths - 1; % Include LOS as one of the paths
%nants = size(bsteer,1);
std = 2; % in degrees
std = std*180/pi; % in radians
phi = rand([1,npaths,nsamples])*2*pi;
phi = repmat(phi,[nants,1,1]);
theta = std*randn([nsamples,npaths]) + losdir;
theta = theta*pi/180;

outsteer = repmat(bsteerdir(losdir,nants),1,nsamples) + squeeze(sum(exp(j*phi) .* bsteerdir(theta,nants),2));



