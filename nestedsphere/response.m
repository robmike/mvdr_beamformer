function s = response(coeffs, center, radius,varargin)
% Analyze the magnitude of steering vectors located in 
% concentric hyperballs with radiuses given in vector 'radius' and center
% 'center' after being filtered by (FIR) filter given by coeffs.

if ~isreal(center) || ~isreal(radius)
   error 'Real arguments required for radius and center' 
end
if any(radius<0)
   error 'Radius must be positive' 
end

bins = [];
normalize = false;

if nargin>3
    if isscalar(varargin{1});
        nummagbins = varargin{1};
    else bins = varargin{1};
    end
    if nargin>4
    if strcmp(varargin{2}(1:length('normalize')),'normalize')
       normalize = true; 
    end
    end
end

nradius = length(radius);
n = length(center);

s = cell(1,nradius);

% % Ratio of volumes of a n-hypercube and hypersphere inscribed in said
% % hypercube is given
% % for i=1:10
% %     v=volhypsphere(radius,i);
% %   ratio = v/(2*radius)^i*100;
% %    fprintf('%i: %.2f %.2f\n',i,v,ratio);
% % end
% % return;
% % ratio = volhypsphere(radius,n)/(2*radius)^n*100
% % ??This provides an estimate to the total number of points that will be in
% % the hypersphere inscribed in the hypercube. Estimate improves as density
% % of points increases.
% % (Does it?)
% npoints = 1e5; % Desired number of points in hypercube
% %npoints = npoints^(1/n); % Number of points per dimension
% npoints = 2^round(log2(npoints)/n)% Nearest power of two: points/dimension
% res = 2*radius/(npoints-1);
% % res = .1;
% numspheres = length(radius);

% We can't use ndgrid for high dimensional spaces because the largest 
% linear index becomes larger than the maximum variable size
% ogrid = cell(n,1);
% 
% % Construct grid with radius centered at zero
% [ogrid{1:n}] = ndgrid(-radius:res:radius);
% ogrid = cell2mat{ogrid}

numsamples = 1e5; % Total number of random samples
blocksize = 1e5; % Number of samples to simulate per iteration
numruns = ceil(numsamples/blocksize);


for k = 1:nradius
for i = 1:numruns
m = i*blocksize;
% Uniformly distributed points in hyperball centered at zero

x = randsphere(m,n,radius(k))';

%Sanity check
if any(vnorm(x,1)>radius(k)+1e-3) % All points should be within hypersphere
   warning('Normx>radius (Points not all within hypersphere)') 
   fprintf('Normx=%.2g > radius=%.2g\n',norm(x),radius(k));
end

% Recenter to given center
x = x + repmat(center,1,size(x,2));

% FIXME: Normalization doesn't distribute points uniformly in hypersphere,
% but without it energy is injected/remove from the system
% This doesn't even generate points uniformly on portion of the surface 
% of the hypersphere with norm sqrt(n) that is contained within the 
% specified hyperball centered at 'center' (at least I don't know that
% it does)
if normalize
    x = x./repmat(vnorm(x,1),n,1); 
end

 % Magnitude of output of filter when excited by input 'x'
if(~isreal(coeffs) && 2*size(coeffs,1)==n)    
 y = ri2complex(x);
else
    y=x;
end

if size(coeffs,1)~=size(y,1)
   error('Size mismatch between coeffs and dimension of space');
end

mag = abs(coeffs'*y);

if isempty(bins)
   magbins = linspace(min(mag),max(mag),nummagbins+1);
else
    magbins= bins;
end

h = histc(mag,magbins);

if i==1
    s{k}.hist = h;
else
    s{k}.hist = h+s{k}.hist;
end
s{k}.magbins = magbins;

end
end

end

%------------------------------------------------------
function v=volhypsphere(r,n)
% Compute volume of n dimensional hypersphere

v=pi^(n/2)*r.^n/gamma(n/2+1);

% Another method
% if rem(n,2)==0 % Case n even
%     v=pi^(n/2)*r^n/factorial(n/2);
% else
%     v=2^( (n+1)/2 )* pi^( (n-1)/2 )*r^n/factd(n);
% end 
end
