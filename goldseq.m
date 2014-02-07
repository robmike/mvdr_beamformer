function gold=goldseq()
% Generate gold sequences of length
% Later: extend this to arbitrary length, preferred sequences... etc.

m=5; % number of bits in LFSR (deg of LFSR polynomial)
n=2^m-1; % length of maximal length sequence


% LFSR coefficients are taken from Matlab help page for the Simulink Gold Sequence generator block (Why is there no equivalent Matlab function?)

% powers of lfsr polynomial coefficients that are non-zero (x^3+1 -> [3 0])
lfsr1 = [5 2 0];
lfsr2 = [5 4 3 2 0];

% Initial state of the LFSRs 
% We should be able to use this as first m bits of the generated preferred sequences but 'filter' does not seem to output these bits 
init1 = [1 1 1 0 0];
init2 = [1 0 0 1 1];

% convert to polynomial representation ([3 0] -> [1 0 0 1] == x^3+1)
gen1 = zeros(1,m+1); % generator polynomial
gen2 = zeros(1,m+1);

gen1(lfsr1 + 1) = 1;
gen1=fliplr(gen1);
gen2(lfsr2 + 1) = 1;
gen2=fliplr(gen2);

% Input nothing ('0') into the filter obtaining zero-input response
% Conceptually: Load initial state and then switch off input, letting LFSR run through feedback alone
% Feedforward coefficients are arbitrary since there is no input
pref1 = filter(1,gen1,zeros(1,n),init1);
pref2 = filter(1,gen2,zeros(1,n),init2);

% Working in GF2 so we need to use mod 2. Can do this at the end because filter is *linear*
% We may be able to instead declare input as elements of GF(2) but this requires comm. toolbox.
pref1 = mod(pref1,2); 
pref2 = mod(pref2,2);

% Generate Gold sequences from preferred sequences
gold=zeros(n+2,n);
for i=1:length(pref2)
gold(i,:) = mod(pref1+[pref2(end-i+1:end) pref2(1:end-i)],2);
end
gold(end-1,:)=pref1;
gold(end,:)=pref2;

gold=1-2.*gold; % Convert to binary {-1,1} instead of {0,1}

%{
% Check cyclic auto correlation properties of preferred sequences (should be two valued
corrvec=zeros(1,n);
x=gold(k,:)';
	for i=1:n
		corrvec(i)=unique(sort(x' * [x(i:end); x(1:i-1)] ) );
	end
	unique(sort(corrvec))
	
end
%}


%{
% Check cyclic cross-correlation properties of gold sequences
% Should be three valued {-1,-t(m),t(m)-2}
% t(m)=2^( (m+1)/2 ) + 1 for odd m 
% t(m)=2^( (m+2)/2 ) + 1 for even m 

% Todo.....
%}


