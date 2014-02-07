% Small test to see how norm of error is related to DOA of error
function test_doaerror()
numants = 10;

b = 90;
adoa = ([-b:b])*pi/180; % Actual doa
pdoa = ([-b:b])*pi/180; % Presumed doa


[x,y] = meshgrid(pdoa,adoa);
[m,n] = size(x);

derr = x-y;
e = pi*(sin(x)-sin(y));

w = zeros(1,1,numants);
w(1,1,:) = [0:numants-1];
%t = zeros(m,n,numants);
t = repmat(w,[m,n,1]);

z = sum(cos(t.*repmat(e,[1,1,numants])),3);
z = sqrt(2*(numants-z));
%nerror = dfunc(e,numants)
pcolor(x,y,z); shading flat; colorbar

% Plot a cross-section
cvals = [0,30,45,90]*pi/180;
for i=1:length(cvals)
[t,idx(i)] = findclosest(adoa, cvals(i),0.1);
end
figure;
plot(x(idx,:)'*180/pi,z(idx,:)');
end

function nerror = dfunc(adoa,pdoa,numants)
e = pi*(sin(adoa)-sin(pdoa));
nerror = sqrt(2)*sqrt(numants - sum(cos([0:numants-1]*e)));
end