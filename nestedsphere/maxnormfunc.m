function y = extremanormfunc(M,phi)

y = extremanumer(M,phi);
end


function y = extremanumer(M)
% Return extrema in one period (starting at zero) of numerator
% cos((M-1)*t/2)*sin(M*t/2)
% Find extrema (not just maxima)
% n=1:2*M-3;
% y = (2n+1)/(M-0.5)*pi
 n=1:2:(4*M-5);
 y = n/(2*M-1)*pi;
 % First maxima is incorrect due to rapid drop off of csc(x/2). Also
 % maximum at 0 is not accounted for yet.
 y(1) = 0; % This fixes above two issues
end
 
function z = zerosnumer(M)
% Return zeros in one period (starting at zero) of numerator
% cos((M-1)*t/2)*sin(M*t/2)
    n = (0:2*M);
    % for even: n<=2*M  for odd: n<2*M-5/2 (n<=2*M-3)
    
    % These need to be multplied by pi to get actual roots
    % Maximum occurs at 0,M,2M,3M,...
    evroots = 2*(1:2*M-1)/M;
    evroots(M) = []; % Not a root, a maximum
    %evroots = 0;
    oddroots = (2*(0:2*M-3)+1)/(M-1);
    z = unique([evroots oddroots]);
    z = z*pi;   
end