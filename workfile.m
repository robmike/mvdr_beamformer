clear;

randstate = 7;
rand('state',randstate); % So that simulation will be deterministic from run to run
randn('state',randstate);

numchanreal=100; % Number of channel realizations
numfiltreal=10; % Number of filter realizations per channel

numusers=6;
numants=5;
delayspread=3; % Not really delayspread. Actually delayspread = realdelayspread+1. Correct this.
noisevar=1;
chanvar = 1/2; % Variance of channel (Assuming variance is defined as 1/2 E{z*z'})

freq = 900e6; % 900 MHz

transdelay = zeros(numusers,1); % Transmission delay of each user (asynchronism)
carphase = zeros(numusers,1); % Carrier phase of each user
% Total carrier phase of each user (including phase changes due to
% transmission delay)
totcarphase = 2*pi*freq*transdelay - carphase ;


% --------- Signature---------------
%gold=goldseq();
%signature=gold(1:numusers,:);

% Trivial signature (no spreading, no CDMA)
% signature = ones(numusers,1);

%signature = 1-2*randbin(1,31); % Single user spreading gain of 31

% Orthogonal hadamard sequences
%signature=hadamard(16);

% Random linearly independent signatures;
signature=zeros(numusers,16);
while rank(signature)<numusers
   signature=2*randbin(numusers,16)-1;
end

signature=signature(1:numusers,:);

spreadgain=size(signature,2);

signature = signature/sqrt(spreadgain); % Normalize signature energy
% -------------------------------------


% Number of data samples to use in estimating correlation matrix
numdatasamplesbase = numants*(spreadgain+delayspread-1);


% ------ Energy ------
energydb=[1 7 8 9.5 10.5 12]; % energy of users in dB
energydb=energydb(1:numusers); % Truncate if using less users
energy=fromdB(energydb);

%energydb=ones(1,numusers);
%energy=fromdB(energydb);


user1snrdb=0:12; % SNR values (in dB) of user of interest to plot
user1snr=fromdB(user1snrdb);
% ---------------------

totruns=length(user1snr)*numchanreal*numfiltreal;


startcpu=cputime;
startwall=clock;

progressbar(0);

samplesmult = [1 2 5 10];

ber=zeros(1,length(user1snr));
ber_array=zeros(length(samplesmult),length(user1snr));
coeffs=zeros(1,numants*(delayspread+spreadgain-1));


for samplesmultidx = 1:length(samplesmult)

smult = samplesmult(samplesmultidx);
numdatasamples = smult * numdatasamplesbase;

for chanidx=1:numchanreal
% ---- Make channel ----
%channel=ones(numusers,delayspread); % No channel 


% Rayleigh channel (zero mean complex Gaussian)
channel = sqrt(chanvar)*...
    complex(randn(numusers,delayspread),randn(numusers,delayspread)).*...
    exp(-j* (2*pi*freq*repmat(1:delayspread,numusers,1) + ...
    repmat(totcarphase,1,delayspread) ) );

% ????????
% Normalize so that total energy is identical for each path 
% (and equal to unity)
% channel = channel./abs(channel)/numants; % ?Is this right?
% More antennas = more received power
% ???????
% ----------------------
%}


% ----- Beamsteering matrix ------
%beamsteer=ones(numusers,numants,delayspread); % No steering


arrayspace = 1/2; % Arrray element spacing in wavelengths
% Angle of arrival of multipaths in radians. 
% Entry (i,j) is angle of jth multipath of ith user
%mpath_aoa = zeros(numusers,delayspread); % Change to desired AOA
mpath_aoa=pi*rand(numusers,delayspread)-pi/2;

% Matrix multipath AOA augmented to 3 dimensions. third dimension is user,
% first dimension is antenna and third is multipath
mpath_aoa_aug = repmat(mpath_aoa,[1,1,numants]); % numsers-by-delayspread-by-numants
mpath_aoa_aug = permute(mpath_aoa_aug,[1 3 2]); % numusers-by-numants-by-delayspread
beamsteer = exp(  j*2*pi*arrayspace*...
    repmat(0:(numants-1),[numusers,1,delayspread]).* sin(mpath_aoa_aug)   ); 

% Change the indexing because we need to index by  mpath_aoa_aug(:,:,1) 
beamsteer = permute(beamsteer,[2 3 1]); % numants-by-delayspread-by-numusers

% --------------------------------

for snridx=1:length(user1snr)

    energy(1)=user1snr(snridx)*noisevar;

    
% ------ If using get_ber_approx
% Estimate interference correlation and add to noise correlation to obtain
% correlation of interference and noise
 intnocorrmat = estintcorrmat(channel,energy,signature,beamsteer,1000);
 intnocorrmat = intnocorrmat + ...
     noisevar*eye(numants*(delayspread+spreadgain-1));
% ---------------------


bsteer1 = beamsteer(:,:,1);
rakevec = rake_coeffs( channel(1,:),energy(1),signature(1,:),...
    bsteer1 );


for filtidx=1:numfiltreal

% Estimate correlation matrix of received data
corrmat = estcorrmat(channel,energy,signature,beamsteer,numdatasamples,noisevar); 
coeffs = mvdr_coeffs(corrmat,rakevec);

%coeffs=rakevec;

while(any(isnan(coeffs)) || any(isinf(coeffs))) % Bad coefficients
    warning('coeferr','Bad coeffs, retrying...\n');
    corrmat = estcorrmat(channel,energy,signature,beamsteer,numdatasamples); 
    coeffs = mvdr_coeffs(corrmat,rakevec);
end


thisber = get_ber(channel, energy, signature, beamsteer,coeffs,noisevar);
%thisber = get_ber_approx(energy(1),coeffs,rakevec,intnocorrmat);
ber(snridx) = ber(snridx) + thisber;


runidx = sub2ind([numfiltreal,length(user1snr),numchanreal],...
    filtidx,snridx,chanidx);

fprintf('MVDR %i of %i: Run %i of %i\n',samplesmultidx,length(samplesmult),runidx,totruns);
pause(0); % Forces fprintf to flush?
progressbar(runidx/totruns );

end % for snridx
end % for chanidx
end % for filtidx

ber=ber/totruns; % Average ber over all runs

ber_array(samplesmultidx,:)=ber;
end % end samplesmult


% Post processing
progressbar(1);

elapsed=cputime-startcpu;
stopwall=clock;


sprintf(...
    'Done.\t Elapsed cpu time: %f minutes\n\t Wall clock time: %f minutes',...
    elapsed/60,etime(stopwall,startwall)/60)

figure;
semilogy(user1snrdb,ber,'o-','markersize',4);
grid on;
