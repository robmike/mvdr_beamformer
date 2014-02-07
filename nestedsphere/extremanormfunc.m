function [extrema,maxima,minima] = extremanormfunc(M,phi)

ubound = floor( (1-sin(phi))*(2*M-1) );
lbound = ceil( (-1-sin(phi))*(2*M-1) );
nu = 3:4:ubound;
nl = fliplr(-3:-4:lbound);
n = [nl nu];
 maxima = n/(2*M-1); 
 nu = 5:4:ubound;
nl = fliplr(-5:-4:lbound);
n = [nl 0 nu];
minima = n/(2*M-1);

maxima = asin( maxima + sin(phi) );
minima = asin( minima + sin(phi) );

% For phi sufficiently far from +-pi/2
% Maximum occurs at sin(theta)-sin(phi) = +- 3/(2*M-1) = x
% and has value sqrt(2*M-1+csc(x*pi/2))
maximum = sqrt(2*M-1+csc(3*pi/2/(2*M-1)));
% First minima occurs at sin(theta)-sin(phi) = +- 5/(2*M-1) = y
% and has value sqrt(2*M-1-csc(y*pi/2))
minimum = sqrt(2*M-1-csc(5*pi/2/(2*M-1)));

% Sometimes, maximum is not achieved on one side.
% We can say (approximately) that maxima is not achieved to the right if
%  sin(phi) > (2*M-4)/(2*M-1)
% Not achieved to the left if
% sin(phi) < (4-2*M)/(2*M-1)
% (Use location of maximum obtained above to show this)

% FIXME: If pi2check is even the maxima at pi/2 and minima at -pi/2.
% If pi2check is odd, minima at pi/2 and maxima at -pi/2
pi2check = M - 1 - ceil((2*M-1)/2*sin(phi));
if maxima(1) < minima(1)
    minima = [-pi/2 minima];
else
    maxima = [-pi/2 maxima];
end
if maxima(end) > minima(end)
    minima = [minima pi/2];
else
    maxima = [maxima pi/2];
end

extrema = sort([maxima minima]);
end


function y = extremanumer(M)
% Return extrema in one period (starting at zero) of numerator
% cos((M-1)*t/2)*sin(M*t/2)
% Find extrema (not just maxima)
% n=1:2*M-3;
% y = (2n+1)/(M-0.5)*pi
 n=1:2:(4*M-5);
 y = n/(2*M-1); % This should be *pi but we will work with integers
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