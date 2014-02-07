clear;

corr_method = 'ideal'; % estimate, ideal exact
coeff_method = 'mvdr'; % auxvec, mvdr, rake

% Force simulation to be pseudo-random but 'repeatable'
randstate = 7;
rand('state',randstate); 
randn('state',randstate);

numchanreal=100; % Number of channel realizations

switch coeff_method
    case 'rake'
        numfiltreal=1;
    otherwise
        numfiltreal=10; % Number of filter realizations per channel
end

numauxvecs = 10; % Number of auxiliary vectors

numusers=6;
numants=1;
delayspread=1; % Not really delayspread. Actually delayspread = realdelayspread+1. Correct this.
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
sigsize = 10;
signature=zeros(numusers,sigsize);
while rank(signature)<numusers
   signature=2*randbin(numusers,sigsize)-1;
end

signature=signature(1:numusers,:);

spreadgain=size(signature,2);

signature = signature/sqrt(spreadgain); % Normalize signature energy
% -------------------------------------



% Number of data samples to use in estimating correlation matrix

switch corr_method
    case 'estimate'
        numdatasamples_array = 76;%[1 2 5 10]*numants*(spreadgain+delayspread-1); %[1 2 5 10]
    case 'ideal'
        numdatasamples_array = 0;
    otherwise
        error('Unknown corrmethod %s',corr_method);
end


if(delayspread==1)
    numbits=1; % Number of bits that affect the bit of interest
else
    numbits=3;
end
% Generate all combinations of information bits of interfering users and 
% infobit of interest of user of interest fixed (to -1)
% Combinations arranged in columns each row repeated numbit times
nbitcombos=2^(numusers*numbits-1);
bitcombos = [-ones(nbitcombos,1) 2*de2bi(0:nbitcombos-1)-1]';
%infobit_array(numusers,numbits);

% ------ Energy ------
%energydb=[1 7 8 9.5 10.5 12]; % energy of users in dB
energydb=20*ones(1,6);
energydb=energydb(1:numusers); % Truncate if using less users
energy=fromdB(energydb);

%energydb=ones(1,numusers);
%energy=fromdB(energydb);


user1snrdb=0:12; % SNR values (in dB) of user of interest to plot
user1snr=fromdB(user1snrdb);
% ---------------------

% Total number of loop iterations. *Not* runs for a given method and SNR
totruns=length(user1snr)*numchanreal*numfiltreal*length(numdatasamples_array);


startcpu=cputime;
startwall=clock;

progressbar(0);


ber_array=zeros(length(numdatasamples_array),length(user1snr));
ber_uncond_array=ber_array;
coeffs=zeros(1,numants*(delayspread+spreadgain-1));


sysparams.numusers = numusers;
sysparams.delayspread = delayspread;
sysparams.spreadgain = spreadgain;
sysparams.numants = numants;
sysparams.noise = noisevar;

for chanidx=1:numchanreal

%channel = raylchan(numusers,delayspread,freq,totcarphase,chanvar);
channel = ones(numusers,delayspread); % No channel

%beamsteer = bsteeruni(numusers,delayspread,numants);
beamsteer = ones([numants,delayspread,numusers]);



for samplesidx = 1:length(numdatasamples_array)

    numdatasamples = numdatasamples_array(samplesidx);


    for snridx=1:length(user1snr)

        energy(1)=user1snr(snridx)*noisevar;
        bsteer1 = beamsteer(:,:,1); % Beamsteering array for user 1
        
        sigarray = get_sigarray(channel,energy,signature,beamsteer);
        rakevec = sigarray(1,:)';
        
        coeffs = rakevec; % If using other filter, these will get overwritten

        for filtidx=1:numfiltreal

            % Estimate correlation matrix of received data
            switch corr_method
                case 'estimate'
                    corrmat = estcorrmat(sigarray,numdatasamples,noisevar,sysparams);
                case {'ideal','exact'}
                    corrmat = idealcorr(sigarray,noisevar);
                otherwise
                    error('Unknown correlation method, %s', corr_method);
            end


            switch coeff_method
                case 'rake'
                    % Do nothing we've already calculated rake coefficients
                case 'mvdr'
                    coeffs = mvdr_coeffs(corrmat,rakevec);
                case 'auxvec'
                    [coeffs,coeffs_uncond]=...
                        auxvec_coeffs(corrmat,rakevec,numauxvecs);
                otherwise
                    error('Unknown coeff function %s', coeff_method);
            end

            % Calculate BER

            if strcmp('auxvec',coeff_method)
                thisber = get_ber(sigarray, bitcombos,coeffs,noisevar);
                ber_uncond_array(samplesidx,snridx) = ber_uncond_array(samplesidx,snridx) + thisber;
            end

            thisber = get_ber(sigarray, bitcombos,coeffs,noisevar);
            ber_array(samplesidx,snridx) = ber_array(samplesidx,snridx) + thisber;

        end % for filtidx

        runidx = sub2ind([numfiltreal,...
            length(user1snr),length(numdatasamples_array),numchanreal],...
            filtidx,snridx,samplesidx,chanidx);

        fprintf('Run %i of %i: %.1f%%\n',runidx,totruns,100*runidx/totruns);
        pause(0); % Lets us control-c more easily
        progressbar(runidx/totruns);

    end % for snridx
end % for samplesidx
end % for chanidx


numtrials = filtidx*chanidx; % Trials at a given SNR (and fixed # of samples)

ber_array=ber_array/numtrials; % Average ber over all runs
if strcmp('auxvec',coeff_method)
    ber_uncond_array=ber_uncond_array/numtrials;
end

% Post processing
progressbar(1);

elapsed=cputime-startcpu;
stopwall=clock;

%save auxvec_est_temp.mat

fprintf(...
    'Done.\t Elapsed cpu time: %f minutes\n\t Wall clock time: %f minutes',...
    elapsed/60,etime(stopwall,startwall)/60);

figure;
semilogy(user1snrdb,ber_array,'o-','markersize',4);
grid on;
