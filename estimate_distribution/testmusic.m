
N = 1e3;
M = 10;
theta = (0.1)*pi;

ses = bsteer(theta,M);

%ses = imag(ses); % CHange me

P = size(ses,2);
s = ses(:,1);

est = zeros(1,N);
numfft = 512;

nvar = 4;
svar = 1;

nsamples = 100;

srep = sqrt(svar)*repmat(s,1,nsamples);

est = zeros(1,N);

for i = 1:N
    x = srep + sqrt(nvar/2)*complex(randn(M,nsamples), ...
                                    randn(M,nsamples));
    R = x*x'/nsamples;
    [est(i),spectrum] = idealmusic(R,P,numfft);
    if mod(i,1000)==0
        fprintf('%i of %i (%.2f%%)\n',i,N,i/N*100);
        %fflush(1);
    end
end


%unique(est)*2
hist(est*2,200);
var(est*2)
%anderson_darling_test(est,'normal')