tic;

clear; close all;
load('c:\autosavesims\20060830T195357-snr.mat');
legtext = {cstr{:},'MVDR, with mismatch','MVDR, no mismatch','MVDR, no mismatch (Ideal $\mathbf{R}$)'};
foosnr;


clear; close all;
load('c:\autosavesims\20060831T005354-samples.mat');
legtext = {cstr{:},'MVDR, with mismatch','MVDR, no mismatch','MVDR, no mismatch (Ideal $\mathbf{R}$)'};
foosamples;


clear; close all;
load('c:\autosavesims\20060830T213243-direrr.mat');
legtext = {cstr{:},'MVDR, with mismatch','MVDR, no mismatch','MVDR, no mismatch (Ideal $\mathbf{R}$)'};
foodirerr;

toc