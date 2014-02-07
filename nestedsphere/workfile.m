M = 10
t=0:0.01:2*pi;
%%
% figure;
% plot(t,[cos( (M-1)/2*t ).* sin( M/2*t );cos( (M-1)/2*t ); sin( M/2*t )])
%%

f = @(x) sin(x/2) + sin((M-1/2)*x); % f has a period 4*pi
phi = 72*pi/180;%asin((2*M-4)/(2*M-1)); %60*pi/180;
[extrema,maxima,minima] = extremanormfunc(M,phi);
%maxima = z(1:2:end); minima = z(2:2:end);
h = @(x) M*cos((M-1)/2*x).*diric(x,M); %sqrt(2*M-2*f(x));
% sin(t/2); -sin( (M-0.5)*t )
%figure;
%plot(t, [h(t)],maxima,h(maxima),'+',minima,h(minima),'x');
g = @(x) cos((M-1)*x/2).*sin( M*x/2 )./sin(x/2);
%plot(t,g(t),maxima,g(maxima),'^',minima,g(minima),'v')

%fg = @(x) M - normfunc(0,x,M).^2/2;
nfunc = @(x) normfunc(phi,x,M);
t = -pi:0.01:pi;
figure;
plot(t,nfunc(t),maxima,nfunc(maxima),'+',minima,nfunc(minima),'x')
hold on
h = line([-pi/2 pi/2;-pi/2 pi/2],[0 0; 5 5])
set(h,'linestyle',':')
hold off;

%% Test tranformation from w to theta
f = @(ed,M)  sqrt(2*M-1-(2*M-1)*diric(ed,(2*M-1)))
x=pi*(-sin(phi)-1):0.01:pi*(-sin(phi)+1);
figure; 
plot(asin(x/pi+sin(phi)),f(x,10));
%%
taui = 0:0.01:asin(3/(2*M-1)+sin(phi));

% Choose the standard deviation as distance to first maximum
eststd = 2*(asin( 3/(2*M-1) + sin(phi) ) - phi);
%figure;
%h = plot(taui,2*taui/pi,taui,erfc(-taui/sqrt(2)/eststd)-1);
%legend('Int','presumed')
%% lkjljklkj

randstate = 7;
rand('state',randstate);
randn('state',randstate);

maxcpfreq = 30*60; % Don't checkpoint more often than this (in seconds)
lastcheckpoint = clock;


chantype = 'none'; % coscat, incscat, none, nearfield, inhomo
enablehist = false;
estimatecorrmat = false; % Set true to simulate SMI and ideal
relativeint = false; % Are interferer DOAs relative to presumed SOI?

logname = 'socp.log';
logfile = fopen(logname,'w');

numants = 10;
nbins = 500; % Number of histogram bins

% min(numsamples_array) must be at least numants
if estimatecorrmat
    numsamples_array = [3]*numants; %[1 2 5 10]*numants;
    numchanreal = 10; %Number of realizations to average over
else
    numsamples_array = [];
    numchanreal = 1;
end

numsamplesets = length(numsamples_array);

soidir_array = [3]*pi/180; % Incident angle of signal of interest.
fprintf('soidir is ');
fprintf('%f ',soidir_array*180/pi); % A check to ensure degrees/rads are not messed up
fprintf('\n');
soisnrdb = -10:20;
soisnr = fromdB(soisnrdb);

numsoidir = length(soidir_array);

if relativeint && numsoidir~=1
   error('You did not implement relative interferers while sweeping soidir yet'); 
end

noisevar = 1; % Noise power

%intdir = [30,50]*pi/180;
% relintdir is interferer direction relative to presumed SOI dir
% if relativeinterferers == true, otherwise absolute interferer directions
% FIXME: Insert a check to make sure absolute DOA is in [0,pi]?
relintdir = [30,50]*pi/180;%[-30,30,40,50,70]*pi/180; % Interferer directions relative to presumed
numint = length(relintdir);
intsnrdb = 30*ones(1,numint); % Per sensor SNR in dB
intsnr = fromdB(intsnrdb);



direrr_array = [-4]*pi/180; % Can be negative also
numdirerr = length(direrr_array);

errnorm = zeros(numsoidir,numdirerr);

% On what index values should we calculate histogram.
%snr_t = linspace(1,length(soisnr),4);
%err_t = linspace(1,numdirerr,4);

% for i = 1:numint
%     temp = bsteerdir(intdir(i),numants);
%     sigs(:,i+1) = temp; %/norm(temp);
% end


intloc = cell(numsoidir,numdirerr);
int2estdist = zeros(numsoidir,numdirerr,numint);
int2actdist = zeros(numsoidir,numdirerr,numint);





% $$$ % Testing to see how far various vectors are from actual steering vector
% $$$ for i = [0:5]*pi/180
% $$$     fprintf('For %.2f degrees int diff norm is %.2f\n',i*180/pi, ...
% $$$             norm(bsteer-bsteerdir(i,numants)) )
% $$$ 
% $$$ end




% intdiffnorm = vnorm(repmat(bsteer,1,length(intdir)) - sigs(:,2:end),1);
% fprintf(['Norm of diff between interferers and beamsteering vector: ' ...
%          '%.2f\n'],intdiffnorm);

% Norm of difference between each interferer and each presumed steering
% vector. Entry (i,j) is the norm of difference between ith bsteer estimate
% and jth interferer
% repsigs = repmat(sigs(:,2:end),[1,1,length(bsteerest_array)]);
% repsigs = permute(repsigs,[1 3 2]); % numants-by-nbsteerest_array-by-nintdir
% intdiffnorm_presumed = vnorm(repmat(bsteerest_array,[1,1,length(intdir)])...
%     - repsigs,1);
% clear repsigs;




% Constraints
con{1}.atten = [1 1e-2];
con{1}.radius = [1 3];

%con{2}.atten = 1; 
%con{2}.radius = con{1}.radius(2);

%con{3}.atten = 1; 
%con{3}.radius = con{1}.radius(1);

%con{4}.atten = [1 0.9 1e-2];
%con{4}.radius = [1 2 3];


numscen = length(con); % Number of scenarios

cstr = cell(1,numscen);
for i = 1:numscen
    cstr{i} = sprintf('>=%.3g in r<%.2f', con{i}.atten(1),con{i}.radius(1));
    for k=2:length(con{i}.atten)
        cstr{i} = [cstr{i} sprintf(', >=%.2f in r<%.2f', con{i}.atten(k), ...
                                   con{i}.radius(k))];
    end
end
fprintf('%s :: ',cstr{:});
fprintf('\n');

oatten = 0.1:0.1:1;


socpsinr = zeros(numscen,length(soisnr),numdirerr,numsamplesets+1,numsoidir);
%socp.sinr.ideal = zeros(numscen,length(soisnr),numdirerr);

mvdrsinr = zeros(length(soisnr),numdirerr,numsamplesets+1,numsoidir);
%mvdr.sinr.ideal = zeros(length(soisnr),numdirerr);

omvdrsinr = zeros(length(soisnr),numdirerr,numsamplesets+1,numsoidir);
%omvdr.sinr.ideal = zeros(1,length(soisnr));

%socpw = zeros(numscen,length(soisnr),numdirerr,numants);
socpw = zeros(numscen,length(soisnr),numdirerr,numsamplesets+1,numants,numsoidir);

%mvdr.w = cell(length(soisnr),numdirerr);
mvdrw = cell(length(soisnr),numdirerr,numsamplesets+1,numsoidir);

%omvdr.w.ideal = cell(1,length(soisnr));
omvdrw = cell(1,length(soisnr),numdirerr,numsamplesets+1,numsoidir);


% FIXME: Fix histogram stuff to deal with simulation of multiple sample
% sizes
nestedsinr_bsteer = zeros(numscen,length(soisnr),numdirerr,length(oatten));
nestedsinr_int = nestedsinr_bsteer;
nestedsinr_out = nestedsinr_bsteer;


% Generate random information bits
randinfobits_array = 2*floor(2*rand([numint+1,max(numsamples_array),numchanreal]))-1;

%----------------- Simulation loop starts here ---------------
starttime = clock;
simdate = datestr(now,30);

for sdiridx = 1:numsoidir

soidir = soidir_array(sdiridx);
tbsteer = bsteerdir(soidir,numants); % 'True' (distortionless), unit energy steering vector
 
direst_array = soidir + direrr_array;
bsteerest_array = bsteerdir(direst_array',numants);

switch chantype
    case 'coscat',
        bsteer = bsteer_coscat(soidir,numants,5,1); % coherent scattering
    case 'none',
        bsteer = bsteerdir(soidir,numants); % Beamsteering vector lossless channel
    case 'incscat',
         bsteer = bsteer_incscat(soidir,numants,5,max(numsamples_array));
    case 'nearfield',
         bsteer = bsteer_nearfield((numants-1)^2/4,numants);
    case 'inhomo',
         bsteer = bsteer_inhomo(soidir,numants);
    otherwise,
        error('Unknown chantype')
end

bsteer_ri = [real(bsteer); imag(bsteer)];
%bsteer = bsteer/norm(bsteer);

switch chantype
    case {'none','coscat','nearfield'}
% Norm of steering vector error
errnorm(sdiridx,:) = vnorm(repmat(bsteer,1,length(direrr_array)) - ...
                         bsteerest_array , 1);
end

for erridx = 1:numdirerr
    intdir = relintdir;
    if relativeint
       intdir = intdir + direst_array(:,erridx);
    end
    intloc{sdiridx,erridx} = bsteerdir(intdir,numants);
%int2estdist(sdiridx,erridx,:) = vnorm(repmat(bsteerest_array(:,erridx),1,length(intdir)) - intloc{sdiridx,erridx},1);
%int2actdist(sdiridx,erridx,:) = vnorm(repmat(bsteer,1,length(intdir)) - intloc{sdiridx,erridx},1);
end


% sigs = [bsteer zeros(numants,numint)];
%sigs = zeros(numants,numint+1);

[t,snridx] = findclosest(soisnr,todb(10))

soivar = soisnr(snridx)*noisevar;


for erridx = 1:numdirerr
for chanidx = 1:numchanreal
    
 
bsteerest = bsteerest_array(:,erridx); 
bsteerest_ri = [real(bsteerest);imag(bsteerest)];

%sigs(:,2:end) = intloc{sdiridx,erridx};
intsigarray = repmat(sqrt(noisevar*intsnr),numants,1).*intloc{sdiridx,erridx};
tsigarray = [sqrt(soivar)*tbsteer intsigarray];
%sigarray = sigs.*repmat([sqrt(soivar) sqrt(noisevar*intsnr)],numants,1);

% Generate random received vectors samples for maximum numsamples
if estimatecorrmat
randinfobits = randinfobits_array(:,:,chanidx);
switch chantype
    case {'none','coscat','nearfield','inhomo'}
        randsamples = [sqrt(soivar)*bsteer intsigarray]*randinfobits;
    case 'incscat'
        randsamples = intsigarray*randinfobits(2:end,:) + sqrt(soivar)*bsteer.*repmat(randinfobits(1,:),numants,1);
    otherwise
        error('unknown channel type')
end
% Add noise
randsamples = randsamples + sqrt(noisevar)*randn(numants,max(numsamples_array));
end

intnocorrmat = idealcorr( intsigarray,noisevar );
    
for sampleidx = 1:numsamplesets+1
    
    if estimatecorrmat && ( sampleidx~=(numsamplesets+1) )
    numsamples = numsamples_array(sampleidx);
    samples = randsamples(:,1:numsamples); % Use appropriate number of samples
    end
    
    if sampleidx==numsamplesets+1
        corrmat = idealcorr(tsigarray,noisevar);
    else
        corrmat = 1/numsamples * samples * samples';
    end



corrmat_real = real(corrmat);
corrmat_imag = imag(corrmat);
R = [corrmat_real -corrmat_imag; corrmat_imag corrmat_real];


for scenidx = 1:numscen 

y = socp_coeffs(bsteerest_ri,con{scenidx}.atten,R,con{scenidx}.radius,logfile);
objfunc = @(x) trashobjfun([1 x],bsteerest_ri,con{1}.atten,con{1}.radius,R,tbsteer,soivar,intnocorrmat)

[optval,fval] = fmincon(objfunc,[0.1 1e-3],[-1 1],[0],[],[],[0;0],[1;1])



coeffs = complex(y(1:end/2),y(end/2+1:end));

if ~abs(imag(coeffs'*bsteerest)<1e-3) % Testing if Vorobyov and
                                        % Lorenz are same
    warning('SOCP:imagzero','Imaginary part of inner product not zero');
end


socpw(scenidx,snridx,erridx,sampleidx,:,sdiridx) = coeffs;
socpsinr(scenidx,snridx,erridx,sampleidx,sdiridx) = get_sinr(coeffs,tbsteer,soivar,intnocorrmat)/chanidx + ...
    socpsinr(scenidx,snridx,erridx,sampleidx,sdiridx)*(chanidx-1)/chanidx;

if enablehist
    fprintf('scenidx=%i,snridx=%i erridx=%i\n',scenidx,snridx,erridx)
        histsphere{scenidx,snridx,erridx} = ...
            response(coeffs,bsteerest_ri,histradiuses,nbins);
        nhistsphere{scenidx,snridx,erridx} = ...
            response(coeffs,bsteerest_ri,histradiuses,nbins,'normalize');
        shistsphere{scenidx,snridx,erridx} = ...
                sphereresponse(coeffs,bsteerest_ri,histradiuses,nbins,nbins);
end


end % for scenidx=1:

w = mvdr_coeffs(corrmat,bsteerest);
mvdrw{snridx,erridx,sampleidx,sdiridx} = w;
mvdrsinr(snridx,erridx,sampleidx,sdiridx) = get_sinr(w,tbsteer,soivar,intnocorrmat)/chanidx + mvdrsinr(snridx,erridx,sampleidx,sdiridx)*(chanidx-1)/chanidx;
if enablehist
    fprintf('snridx=%i erridx=%i\n',snridx,erridx)
    histmvdr{snridx,erridx} = response(w,bsteerest_ri,histradiuses,magbins);
    nhistmvdr{snridx,erridx} = ...
            response(w,bsteerest_ri,histradiuses,magbins,'normalize');
    shistmvdr{snridx,erridx} = ...
                sphereresponse(w,bsteerest_ri,histradiuses,nbins,nbins);
end

% Note that since the interferer directions depend on presumed steering
% vector (i.e. erridx) then so does ideal MVDR filter coefficients
switch chantype,
    case 'incscat',
        % Principle eigenvector method
        Rs = tbsteer*tbsteer';
        [U,D] = svd(soivar*intnocorrmat\Rs);
        w = U(:,1); % Principle eigenvector
        omvdrw{snridx,erridx,sampleidx,sdiridx} = w;
        omvdrsinr(snridx,erridx,sampleidx,sdiridx) = real(w'*Rs*w)/real(w'*intnocorrmat*w);
    case {'none','coscat','nearfield','inhomo'}
        w = mvdr_coeffs(corrmat,bsteer);
        omvdrw{snridx,erridx,sampleidx,sdiridx} = w;
        omvdrsinr(snridx,erridx,sampleidx,sdiridx) = soivar*real(tbsteer'*inv(intnocorrmat)*tbsteer);
    otherwise,
        error('unknown chantype')
end

end % for samplesidx=1:

if estimatecorrmat
    progresstext([chanidx erridx snridx sdiridx], [numchanreal numdirerr length(soisnr) numsoidir]);
end

end % for chanidx = 1:


% ---- Save a checkpoint
if etime(clock,lastcheckpoint) > maxcpfreq
lastcheckpoint = clock;
save('c:\autosavesims\checkpoint.mat','-v6');
fprintf('Saved checkpoint at %s\n',datestr(now,13));
end

end % for erridxidx=1:

if ~estimatecorrmat
    progresstext([snridx sdiridx], [length(soisnr) numsoidir]);
end




if enablehist
    fprintf('snridx=%i erridx=%i\n',snridx,erridx)
    histopt{snridx} = response(w,bsteer_ri,histradiuses,magbins);
    nhistopt{snridx} = ...
            response(w,bsteerest_ri,histradiuses,magbins,'normalize');
    shistopt{snridx} = ...
                sphereresponse(w,bsteerest_ri,histradiuses,nbins,nbins);
end


end % for sdiridx

fprintf('Creating plots...\n');
fprintf('Done. Total wall time: %.2f min\n', etime(clock,starttime)/60);
datefin = datestr(now);
fprintf('Finished at: %s\n',datefin);

fclose(logfile);

% FIXME: Create autosavesims directory if it doesn't exist
 save(['C:\autosavesims\' simdate '.mat'])
 fprintf('Auto-saved workspace to %s\n',[simdate '.mat'])

 nelliplot


