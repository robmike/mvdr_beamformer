function [c11 C] = MUSIC_variance(w,numants)

numsamples = 30;
nvar = 4;
svar = 1;
numusers = length(w);

N = numsamples;
w = torow(w);

A = sqrt(svar/numants)*bsteerdirfreq(w,numants);
D = sqrt(svar/numants)*diffbeamsteer(w,numants);

P = eye(numusers);
Pinv = inv(P);
H = D'*(eye(numants) - A*(A'*A)^(-1)*A')*D;
hid = inv(H.*eye(numusers));
C = real(H.* (Pinv + sqrt(nvar)*Pinv*inv(A'*A)*Pinv).');
C = sqrt(nvar)/(2*N)*hid*C*hid;

c11 = C(1,1);

end

function z = bsteerdirfreq(w,numants)
x = j*(0:numants-1)';
x = repmat(x,1,length(w));
w = repmat(w,numants,1);
z = exp(x.*w);
end

function dz = diffbeamsteer(w,numants)
% Derivative of beamsteering vector with respect to electrical steering
% angle w. w may be a vector, in which cas 

x = j*(0:numants-1)'; 
w = torow(w);
x = repmat(x,1,length(w));
w = repmat(w,numants,1);

dz = x.*exp(x.*w);

end


