function bsteer = bsteer_nearfield(distance,numants)
% Distance in wavelengths
% Source placed at center of array 
% Vorobyov: distance of (M-1)^2/4 wavelengths

M = numants;
i = (0:M-1)';
%bsteer = exp(j*pi*sqrt( (M-1)^4/4 + ((M-1)/2 - i).^2))'
bsteer = exp(-j*2*pi*sqrt( distance.^2 + ((M-1)/2 - i).^2));

end