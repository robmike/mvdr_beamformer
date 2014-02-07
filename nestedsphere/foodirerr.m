% Script to plot data obtained after running nestedellipse.m
close all

plotdir = 'plotsbb/';
lfontsize = 14; % Font size for legend
xfontsize = 14;
yfontsize = 14;

texinterp = 'none'; % Tex interpreter

epsexport = false;

% Create plot direcotory if it does not exist
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end

ytext = 'Output SINR (dB)';

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

%% Output SINR vs. error in degrees for fixed input SNR
%----------------------------------------
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);
degerror = direrr_array*180/pi;


stitle = sprintf('SNR=%.2f dB',snrx);
hf = plot(degerror,todB(squeeze(socpsinr(:,sidx,:,ridx))),'-+',...
        degerror, todB(mvdrsinr(sidx,:,ridx)),'--x',...
        degerror, todB(repmat(mvdrsinr(sidx,zidx,ridx),1,length(degerror))),...
        degerror,todB(omvdrsinr(sidx,:,ridx)),'--');

cyclelines(gca);
%set(hf(5),'visible','off'); % Hide MVDR no mismatch plot
delete(hf(5))
% grid on;
%title(stitle);
ylabel(ytext,'fontsize',yfontsize);
xlabel('Steering error in degrees','fontsize',xfontsize);
h=legend(legtext{[1:4 6]},'Location','Southeast');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize,'interpreter',texinterp);



% Increase y axis top margin
c = 0.08;  % 8 percent
% Computes position of new axis limit so that fraction of 'new part of
% axis' over total size of the new axis size is c
logfrac = @(x,y,c) y^(1./(1-c))*x^(-c./(1-c));
yl = ylim;
yl(2) = logfrac(yl(1),yl(2),c);
ylim(yl);

fstr = sprintf('error-sinr%.0f',snrx);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir '/' fstr '.eps'])
end

polishfigs
close all;

plotdir = 'plotsla/';
texinterp = 'latex';
%% Output SINR vs. error in degrees for fixed input SNR
%----------------------------------------
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);
degerror = direrr_array*180/pi;


stitle = sprintf('SNR=%.2f dB',snrx);
hf = plot(degerror,todB(squeeze(socpsinr(:,sidx,:,ridx))),'-+',...
        degerror, todB(mvdrsinr(sidx,:,ridx)),'--x',...
        degerror, todB(repmat(mvdrsinr(sidx,zidx,ridx),1,length(degerror))),...
        degerror,todB(omvdrsinr(sidx,:,ridx)),'--');

cyclelines(gca);
%set(hf(5),'visible','off'); % Hide MVDR no mismatch plot
delete(hf(5))
% grid on;
%title(stitle);
ylabel(ytext,'fontsize',yfontsize);
xlabel('Steering error in degrees','fontsize',xfontsize);
h=legend(legtext{[1:4 6]},'Location','Southeast');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize,'interpreter',texinterp);



% Increase y axis top margin
c = 0.08;  % 8 percent
% Computes position of new axis limit so that fraction of 'new part of
% axis' over total size of the new axis size is c
logfrac = @(x,y,c) y^(1./(1-c))*x^(-c./(1-c));
yl = ylim;
yl(2) = logfrac(yl(1),yl(2),c);
ylim(yl);

fstr = sprintf('error-sinr%.0f',snrx);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir '/' fstr '.eps'])
end


polishfigs
