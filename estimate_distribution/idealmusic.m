%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: music.m
%
% AUTHOR: Steve Kogon
% 
% DATE: January 18, 1999
%
% DESCRIPTION: This file forms an estimate of the frequency spectrum using 
%              MUSIC algorithm (spectral version), (Schmidt 1986).
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is modified from its original version to accept
% correlation matrix directly
%
%-----------------------------------------------------------
% Copyright 2000, by Dimitris G. Manolakis, Vinay K. Ingle,
% and Stephen M. Kogon.  For use with the book
% "Statistical and Adaptive Signal Processing"
% McGraw-Hill Higher Education.
%-----------------------------------------------------------


function [fest,Rbar] = idealmusic(R,P,Nfft)
% R = ideal correlation matrix
% P = number of complex exponentials
% M = time-window length
% Nfft = number of FFT frequencies

M = size(R,1);
% Compute eigendecomposition and order by descending eigenvalues
[Q0,D] = eig(R);
[lambda,index] = sort(abs(diag(D)),'descend');
Q=Q0(:,index);

% Compute pseudo-spectrum
f = (-Nfft/2:(Nfft/2-1))/Nfft;                % FFT frequencies
Qbar = zeros(Nfft,1);

for n = 1:(M-P)
    Qbar = Qbar + abs(fftshift(fft(Q(:,M-(n-1)),Nfft))).^2;
    %Qbar = Qbar + abs(fftshift(imag(fft(Q(:,M-(n-1)),Nfft)))).^2; %Change
end
Rbar = 1./Qbar;

% Find local maxima (values of Rbar that are larger than their neighbors)
z1 = Rbar(2:(Nfft-1)) - Rbar(1:(Nfft-2));
z2 = Rbar(2:(Nfft-1)) - Rbar(3:Nfft);
peak_index = find((z1 > 0) & (z2 > 0)) + 1;
if Rbar(1) > Rbar(2)
   peak_index = [1 
      peak_index];
end
if Rbar(Nfft) > Rbar(Nfft-1)
   peak_index = [peak_index
                 Nfft];
end
Npeaks = length(peak_index);
f_peaks = f(peak_index);

% Determine the P largest peaks (taken from local maxima)
[dummy,fest_index] = sort(Rbar(peak_index));
fest = f_peaks(fest_index(Npeaks:-1:(Npeaks-P+1))).';
