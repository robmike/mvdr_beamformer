function sigarray = get_sigarray(channel,energy,signature,beamsteer)

spreadgain=size(signature,2); % Number of samples in signature waveform is spreading gain
numusers=size(channel,1);
numants=size(beamsteer,1);
delayspread=size(channel,2); % Largest multipath delay (in samples)

if(delayspread==1)
    numbits=1; % Number of bits that affect the bit of interest
else
    numbits=3;
end

idxboi = (1+numbits)/2;


zeropad=zeros(numants,spreadgain);
sigarray = zeros(numbits*numusers,numants*(spreadgain+delayspread-1));

for i=1:numusers
    sidx = numbits*(i-1)+1;
    % Signature vector of user i after is passed through the channel
    bsteer = beamsteer(:,:,i);
    sout=sigout(channel(i,:),energy(i),signature(i,:),bsteer);
    
    sigarray(sidx,:) = sout(:)';
    
    if(numbits~=1) % Take into account multipath
        
        %Output due to bit i-1 that interferes with bit i:last delayspread-1
        %samples
        temp = [sout(:,end-delayspread+2:end) zeropad]; 
        sigarray(sidx+1,:) = temp(:)';
        
        %Output due to bit i+1 that interferes with bit i:first delayspread-1
        %samples
        temp = [zeropad sout(:,1:delayspread-1)];
        sigarray(sidx+2,:) = temp(:)';
    end
end