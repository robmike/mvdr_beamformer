% Script to plot data obtained after running nestedellipse.m
close all

plotdir = 'plotsbb/';
lfontsize = 14; % Font size for legend
xfontsize = 14;
yfontsize = 14;

ytext = 'Output SINR (dB)';

texinterp = 'none'; % Tex interpreter

epsexport = false;

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

%% Output SINR vs. Input SNR
%----------------------------------------
numplots = 10;
%[nsample,ridx] = findclosest(numsamples_array,3*numants)
[nsample,ridx] = findclosest(numsamples_array,30); % What sample size to use.
ridx = max(ridx,1)

for idx = 1:ceil(numdirerr/numplots):numdirerr;
f = figure; 
degerror = direrr_array(idx)*180/pi;
% stitle = sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
%                  degerror,errnorm(idx));
% stitle = [stitle sprintf(' %.2f',int2actdist(idx,:))];
% stitle = [stitle sprintf(' || intestdist =')];
% stitle = [stitle sprintf(' %.2f',int2estdist(idx,:))];
plot(soisnrdb,todB(socpsinr(:,:,idx,ridx)),'-+',...
    soisnrdb,todB(mvdrsinr(:,idx,ridx)),'--x',...
    soisnrdb,todB(mvdrsinr(:,zidx,ridx)),...
    soisnrdb,todB(omvdrsinr(:,idx,ridx)),'-.'...
    );
%grid on;
%title(stitle);

h=legend(legtext{:},'Location','NorthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR (mismatch)','MVDR (no mismatch).'Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize,'interpreter',texinterp);
xlabel('Input SNR per antenna (dB)','fontsize',xfontsize)
ylabel(ytext,'fontsize',yfontsize)
cyclelines(gca);

fstr = sprintf('snr-sinr%.0f',degerror);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end
end

polishfigs
close all;

plotdir = 'plotsla/';
lfontsize = 14; % Font size for legend
xfontsize = 14;
yfontsize = 14;

texinterp = 'latex'; % Tex interpreter

epsexport = false;

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


[t,zidx] = findclosest(direrr_array,0,1e-6); % Index for no error

%% Output SINR vs. Input SNR
%----------------------------------------
numplots = 10;
%[nsample,ridx] = findclosest(numsamples_array,3*numants)
[nsample,ridx] = findclosest(numsamples_array,30); % What sample size to use.
ridx = max(ridx,1)

for idx = 1:ceil(numdirerr/numplots):numdirerr;
f = figure; 
degerror = direrr_array(idx)*180/pi;
% stitle = sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
%                  degerror,errnorm(idx));
% stitle = [stitle sprintf(' %.2f',int2actdist(idx,:))];
% stitle = [stitle sprintf(' || intestdist =')];
% stitle = [stitle sprintf(' %.2f',int2estdist(idx,:))];
plot(soisnrdb,todB(socpsinr(:,:,idx,ridx)),'-+',...
    soisnrdb,todB(mvdrsinr(:,idx,ridx)),'--x',...
    soisnrdb,todB(mvdrsinr(:,zidx,ridx)),...
    soisnrdb,todB(omvdrsinr(:,idx,ridx)),'-.'...
    );
%grid on;
%title(stitle);

h=legend(legtext{:},'Location','NorthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR (mismatch)','MVDR (no mismatch).'Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize,'interpreter',texinterp);
xlabel('Input SNR per antenna (dB)','fontsize',xfontsize)
ylabel(ytext,'fontsize',yfontsize)
cyclelines(gca);

fstr = sprintf('snr-sinr%.0f',degerror);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end
end

polishfigs;