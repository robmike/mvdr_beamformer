M = 4;
phi = pi/4;
nvar = .264;
Q = @(x) 0.5*erfc(x/sqrt(2));

t = 0:0.01:(2*M-1-csc(5/(2*M-1)*pi/2));

options = optimset('display','off');
a = zeros(1,length(t));
for i=1:length(t)
     a(i) = fminbnd(@(x) (normfunc2(x,0,M).^2 - t(i))^2,0,3*pi/M,options);
end

%a = 0:0.01:pi*(sin(phi)+1)+0.5;

soiCDF = Q(-a/sqrt(nvar)) - Q(a/sqrt(nvar));
%intCDF = sinuniCDF(diff(sin(phi) + [a/pi;-a/pi],1,1));
intCDF = sinuniCDF(sin(phi) + a/pi) - sinuniCDF(sin(phi) - a/pi);

figure;
plot(a,soiCDF,a,intCDF);
legend('soi','int')