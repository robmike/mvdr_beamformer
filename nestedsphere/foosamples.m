% Script to plot data obtained after running nestedellipse.m
close all

plotdir = 'plotsbb/';
lfontsize = 14; % Font size for legend
xfontsize = 14;
yfontsize = 14;

ytext = 'Output SINR (dB)';

texinterp = 'none'; % Tex interpreter

epsexport = false;

% Increase y axis top margin
c = 0.08;  % 8 percent
% Computes position of new axis limit so that fraction of 'new part of
% axis' over total size of the new axis size is c
logfrac = @(x,y,c) y^(1./(1-c))*x^(-c./(1-c));

% Create plot direcotory if it does not exist
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end

% --- Legend strings
cstr = cell(1,numscen);
for i = 1:numscen
    if(length(con{i}.atten)==1)
    cstr{i} = sprintf('1-constraint ($a$=%i)', con{i}.radius(1));
    else
     cstr{i} = sprintf('2-constraint')
    end 
end
fprintf('%s :: ',cstr{:});
fprintf('\n');
legtext = {cstr{:},'MVDR, with mismatch','MVDR, no mismatch','MVDR, no mismatch (Ideal $\mathbf{R}$)'};


[t,zidx] = findclosest(direrr_array,0,1e-6); % Index for no error

[nsample,ridx] = findclosest(numsamples_array,30); % What sample size to use.
ridx = max(ridx,1)

%% Output SINR vs. sample size for fixed input SNR and fixed error
% -----------------------------------------------
if estimatecorrmat
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);

nsa = numsamples_array/numants;

% Find degerror closest to 
[degerror,eidx] = findclosest(direrr_array*180/pi, 3);

% stitle = sprintf('Input SNR=%.1f dB \n',snrx)
% stitle = [stitle sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
%                  degerror,errnorm(eidx))];
% stitle = [stitle sprintf(' %.2f',int2actdist(eidx,:))];
% stitle = [stitle sprintf(' || intestdist =')];
% stitle = [stitle sprintf(' %.2f',int2estdist(eidx,:))];
plot(nsa,todB(squeeze(socpsinr(:,sidx,eidx,1:end-1))),'-+',...
        nsa,todB(squeeze(mvdrsinr(sidx,eidx,1:end-1))),'--x',...
        nsa,todB(squeeze(mvdrsinr(sidx,zidx,1:end-1))),...
        nsa,todB(squeeze(omvdrsinr(sidx,eidx,1:end-1))),'-.')
%grid on;
%title(stitle);



h=legend(legtext{:},'Location','SouthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize,'interpreter',texinterp);
xlabel('Sample size per antenna (samples/antenna)','fontsize',xfontsize)
ylabel(ytext,'fontsize',yfontsize)
cyclelines(gca);

% Increase y axis top margin
c = 0.08;  % 8 percent
yl = ylim;
yl(2) = logfrac(yl(1),yl(2),c);
ylim(yl);

fstr = sprintf('samples-sinr%.0f');
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end

polishfigs;
close all;
end

plotdir = 'plotsla/';
lfontsize = 14; % Font size for legend
xfontsize = 14;
yfontsize = 14;

texinterp = 'latex'; % Tex interpreter
%% Output SINR vs. sample size for fixed input SNR and fixed error
% -----------------------------------------------
if estimatecorrmat
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);

nsa = numsamples_array/numants;

% Find degerror closest to 
[degerror,eidx] = findclosest(direrr_array*180/pi, 3);

% stitle = sprintf('Input SNR=%.1f dB \n',snrx)
% stitle = [stitle sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
%                  degerror,errnorm(eidx))];
% stitle = [stitle sprintf(' %.2f',int2actdist(eidx,:))];
% stitle = [stitle sprintf(' || intestdist =')];
% stitle = [stitle sprintf(' %.2f',int2estdist(eidx,:))];
plot(nsa,todB(squeeze(socpsinr(:,sidx,eidx,1:end-1))),'-+',...
        nsa,todB(squeeze(mvdrsinr(sidx,eidx,1:end-1))),'--x',...
        nsa,todB(squeeze(mvdrsinr(sidx,zidx,1:end-1))),...
        nsa,todB(squeeze(omvdrsinr(sidx,eidx,1:end-1))),'-.')
%grid on;
%title(stitle);



h=legend(legtext{:},'Location','SouthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize,'interpreter',texinterp);
xlabel('Sample size per antenna (samples/antenna)','fontsize',xfontsize)
ylabel(ytext,'fontsize',yfontsize)
cyclelines(gca);

% Increase y axis top margin
c = 0.08;  % 8 percent
yl = ylim;
yl(2) = logfrac(yl(1),yl(2),c);
ylim(yl);

fstr = sprintf('samples-sinr%.0f');
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end

polishfigs;
close all;

end