% Test the extension of the worst case performance optimization
% as presented in Boyd -- Robust Minimum Variance Beamforming
% to the case of multiple nested ellipses, each with different
% constraints on the degree of acceptable attenuation for
% steering vectors that fall within the given ellipsoid.
% Idea: Actual vectors are more likely to be near the
% estimated/presumed steering vector so 'protect' these more
% while protecting steering vector errors further away
% somewhat less.

warning('Not clearing workspace')
len = @(x) length(x); % Python syntax for length command

randstate = 7;
rand('state',randstate);
randn('state',randstate);

warning('off','optim:fmincon:SwitchingToMediumScale');

maxcpfreq = 30*60; % Don't checkpoint more often than this (in seconds)
lastcheckpoint = clock;

chantype = 'none'; % coscat, incscat, none, nearfield, inhomo
findopt = false; % Try and find optimum filter for each SNR point
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
soisnrdb = 10;
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



direrr_array = [-2]*pi/180; % Can be negative also
numdirerr = length(direrr_array);

errnorm = zeros(numsoidir,numdirerr);


intloc = cell(numsoidir,numdirerr);
int2estdist = zeros(numsoidir,numdirerr,numint);
int2actdist = zeros(numsoidir,numdirerr,numint);




% Constraints
con{1}.atten = 1;
con{1}.radius = 3;


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

arange = 0.1:0.05:(sqrt(numants)-0.1);
crange = 10.^(-6:0.005:0); %0:0.005:1; %
[radvals,cvals] = meshgrid(arange,crange);

oneconsinr = 0;
twoconsinr = zeros([size(radvals)]);


optparams = cell(numscen,length(soisnr),numdirerr,numsamplesets+1,numsoidir,2);


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


for snridx = 1:length(soisnr)

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


    % ------------ One constraint --------------------
    y = socp_coeffs(bsteerest_ri,con{1}.atten,R,con{1}.radius,logfile);
    coeffs = complex(y(1:end/2),y(end/2+1:end));
    oneconsinr = get_sinr(coeffs,tbsteer,soivar,intnocorrmat);
    coeffs = complex(y(1:end/2),y(end/2+1:end));
if ~abs(imag(coeffs'*bsteerest)<1e-3) % Testing if Vorobyov and
                                        % Lorenz are same
    warning('SOCP:imagzero','Imaginary part of inner product not zero');
end

% --------- Two constraint --------------
progressbar(0);
for paramidx = 1:length(radvals(:)) 
    y = socp_coeffs(bsteerest_ri,[con{1}.atten cvals(paramidx)],R,[radvals(paramidx) con{1}.radius],logfile);
    coeffs = complex(y(1:end/2),y(end/2+1:end));
    twoconsinr(paramidx) = get_sinr(coeffs,tbsteer,soivar,intnocorrmat);
    coeffs = complex(y(1:end/2),y(end/2+1:end));
if ~abs(imag(coeffs'*bsteerest)<1e-3) % Testing if Vorobyov and
                                        % Lorenz are same
    warning('SOCP:imagzero','Imaginary part of inner product not zero');
end
progressbar(paramidx/length(radvals(:)));
end % paramidx



end % for samplesidx=1:

end % for chanidx = 1:


% ---- Save a checkpoint
if etime(clock,lastcheckpoint) > maxcpfreq
lastcheckpoint = clock;
save('c:\autosavesims\checkpoint.mat','-v6');
fprintf('Saved checkpoint at %s\n',datestr(now,13));
end

end % for erridxidx=1:

% if ~estimatecorrmat
%     progresstext([snridx sdiridx], [length(soisnr) numsoidir]);
% end




if enablehist
    fprintf('snridx=%i erridx=%i\n',snridx,erridx)
    histopt{snridx} = response(w,bsteer_ri,histradiuses,magbins);
    nhistopt{snridx} = ...
            response(w,bsteerest_ri,histradiuses,magbins,'normalize');
    shistopt{snridx} = ...
                sphereresponse(w,bsteerest_ri,histradiuses,nbins,nbins);
end


end % for snridx
end % for sdiridx

fprintf('Creating plots...\n');
fprintf('Done. Total wall time: %.2f min\n', etime(clock,starttime)/60);
datefin = datestr(now);
fprintf('Finished at: %s\n',datefin);

fclose(logfile);

% FIXME: Create autosavesims directory if it doesn't exist
 save(['C:\autosavesims\' simdate '.mat'])
 fprintf('Auto-saved workspace to %s\n',[simdate '.mat'])


% nelliplot
sinr_ratio = twoconsinr/oneconsinr;
figure;
pcolor(radvals,cvals,sinr_ratio); colorbar; shading flat;
figure;
z = sinr_ratio;
z(z>1.1) = 255;
z(abs(z-1)<=0.1) = 127;
z(z<0.9) = 0;
colormap gray;
pcolor(radvals,cvals,z); colorbar; shading flat
caxis([0 255]);
figure; % Logarthmic scale
pcolor(radvals,log10(cvals),sinr_ratio); colorbar; shading flat