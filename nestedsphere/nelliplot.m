% Script to plot data obtained after running nestedellipse.m
plotdir = 'plots/';
lfontsize = 12; % Font size for legend
xfontsize = 12;
yfontsize = 12;

epsexport = false;

% Create plot direcotory if it does not exist
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end

%% Output SINR vs. Input SNR
%----------------------------------------
numplots = 10;
[t,sdiridx] = findclosest(soidir_array,5*pi/180);
%[nsample,ridx] = findclosest(numsamples_array,3*numants)
[nsample,ridx] = findclosest(numsamples_array,30); % What sample size to use.
if isempty(ridx)
   ridx = 1; 
end
ridx = max(ridx,1);

for idx = 1:ceil(numdirerr/numplots):numdirerr;
f = figure; 
degerror = direrr_array(idx)*180/pi;
stitle = sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
                 degerror,errnorm(sdiridx,idx));
stitle = [stitle sprintf(' %.2f',int2actdist(sdiridx,idx,:))];
stitle = [stitle sprintf(' || intestdist =')];
stitle = [stitle sprintf(' %.2f',int2estdist(sdiridx,idx,:))];
semilogy(soisnrdb,socpsinr(:,:,idx,ridx,sdiridx),'-+',soisnrdb,mvdrsinr(:,idx,ridx,sdiridx),'--x',soisnrdb,omvdrsinr(:,idx,ridx,sdiridx),'-.')
%grid on;
title(stitle);

h=legend(cstr{:},'MVDR w/mismatch','Ideal MVDR','Location','NorthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize);
xlabel('Input SNR per antenna (dB)','fontsize',xfontsize)
ylabel('Output SINR','fontsize',yfontsize)
cyclelines(gca);

fstr = sprintf('snr-sinr%.0f',degerror);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end
end



%% Output SINR vs. error in degrees for fixed input SNR
%----------------------------------------
if numdirerr>1
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);
degerror = direrr_array*180/pi;


stitle = sprintf('SNR=%.2f dB',snrx);
semilogy(degerror,squeeze(socpsinr(:,sidx,:,ridx,sdiridx)),'-+',degerror, ...
         mvdrsinr(sidx,:,ridx,sdiridx),'--x', degerror,omvdrsinr(sidx,:,ridx,sdiridx),'--')
% grid on;
title(stitle);
ylabel('Output SINR','fontsize',yfontsize);
xlabel('Steering error in degrees','fontsize',xfontsize);
h=legend(cstr{:},'MVDR w/steering error','Ideal MVDR','Location','Southwest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize);
cyclelines(gca);

fstr = sprintf('error-sinr%.0f',snrx);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir '/' fstr '.eps'])
end
end
%% Output SINR vs. norm of error for fixed input SNR
%----------------------------------------
if numdirerr>1
f=figure;
[nerror_sorted perm] = sort(errnorm(sdiridx,:));

stitle = sprintf('SNR=%.2f dB',snrx);
semilogy(nerror_sorted,squeeze(socpsinr(:,sidx,perm,ridx,sdiridx)),'-+',nerror_sorted, ...
         mvdrsinr(sidx,perm,ridx,sdiridx),'--x', nerror_sorted,omvdrsinr(sidx,perm,ridx,sdiridx),'--')
%grid on;
title(stitle);
ylabel('Output SINR');
xlabel('Norm of steering error');
h=legend(cstr{:},'MVDR w/steering error','Ideal MVDR','Location','Southwest');
set(h,'fontsize',lfontsize);
cyclelines(gca);

fstr = sprintf('nerror-sinr%.0f',snrx);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir '/' fstr '.eps'])
end
end

%% Output SINR vs. sample size for fixed input SNR and fixed error
% -----------------------------------------------
if numsamplesets>1
if estimatecorrmat
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);

nsa = numsamples_array/numants;

% Find degerror closest to 
[degerror,eidx] = findclosest(direrr_array*180/pi, 3);

stitle = sprintf('Input SNR=%.1f dB \n',snrx);
stitle = [stitle sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
                 degerror,errnorm(sdiridx,eidx))];
stitle = [stitle sprintf(' %.2f',int2actdist(sdiridx,eidx,:))];
stitle = [stitle sprintf(' || intestdist =')];
stitle = [stitle sprintf(' %.2f',int2estdist(sdiridx,eidx,:))];
semilogy(nsa,squeeze(socpsinr(:,sidx,eidx,1:end-1,sdiridx)),'-+',nsa,squeeze(mvdrsinr(sidx,eidx,1:end-1,sdiridx)),'--x',nsa,squeeze(omvdrsinr(sidx,eidx,1:end-1,sdiridx)),'-.')
%grid on;
title(stitle);

h=legend(cstr{:},'MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize);
xlabel('Sample size per antenna (samples/antenna)','fontsize',xfontsize)
ylabel('Output SINR','fontsize',yfontsize)
cyclelines(gca);

fstr = sprintf('samples-sinr%.0f',snrx);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end
end
end
%% Output SINR versus SOI DOA for a fixed input SNR and fixed DOA error
if numsoidir>1
f=figure;
[snrx,sidx] = findclosest(soisnrdb,10);
[degerror,eidx] = findclosest(direrr_array*180/pi, 3);
sda = soidir_array*180/pi;
if estimatecorrmat
   ns = numsamples_array(ridx);
else
    ns = inf;
end

stitle = sprintf('Input SNR=%.1f dB \n',snrx);
stitle = [stitle sprintf('%i samples \n',ns)];
stitle = [stitle sprintf('Error = %.2f deg (norm = %.2f)\nintactdist =', ...
                 degerror,errnorm(sdiridx,eidx))];
stitle = [stitle sprintf(' %.2f',int2actdist(sdiridx,eidx,:))];
stitle = [stitle sprintf(' || intestdist =')];
stitle = [stitle sprintf(' %.2f',int2estdist(sdiridx,eidx,:))];
semilogy(sda,squeeze(socpsinr(:,sidx,eidx,ridx,:)),'-+',sda,squeeze(mvdrsinr(sidx,eidx,ridx,:)),'--x',sda,squeeze(omvdrsinr(sidx,eidx,ridx,:)),'-.')
%grid on;
title(stitle);

h=legend(cstr{:},'MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
%h = legend('1 constraint (a=3)','2 constraints','1 constraint (a=1)','MVDR w/mismatch','Ideal MVDR','Location','SouthWest');
set(h,'fontsize',lfontsize);
xlabel('SOI DOA (degrees)','fontsize',xfontsize)
ylabel('Output SINR','fontsize',yfontsize)
cyclelines(gca);

fstr = sprintf('soidoa-sinr%.0f',snrx);
set(f,'userdata',fstr); % Save plot name in Figure.
if epsexport
print('-depsc2',[plotdir fstr '.eps'])
end

end

%% Responses versus outer attenuation
%------------------------------------------
% snridx = find(soisnrdb>=5,1);
% scenidx = 2;
% erridx = 2;
% 
% w = squeeze(nestedsinr_int(scenidx,snridx,erridx,:));
% figure;
% plot(oatten,w,oatten,oatten);
% %grid on;
% 
% 
% hold on;
% wnestnorm = zeros(1,length(oatten));
% crosscorr = zeros(1,length(oatten));
% for i=1:length(oatten)
%     w = wnested{scenidx,snridx,erridx,i};
%     crosscorr(i) = w'*sigs(:,2)/norm(w)/norm(sigs(:,2));
%     wnestnorm(i) = norm(wnested{scenidx,snridx,erridx,i});
% end
% plot(oatten,wnestnorm);
% 
% 
% 
% 
% 
% w = squeeze(nestedsinr_bsteer(scenidx,snridx,erridx,:));
% %figure;
% plot(oatten,w);
% %grid on;
% 
% plot(oatten,abs(crosscorr));
% hold off
% cyclelines(gcf);
% legend('Interference steering','Min attenuation','||w||','Actual steering')
% 
% %%
% w = squeeze(nestedsinr_out(scenidx,snridx,erridx,:));
% figure;
% semilogy(oatten,w);
% grid on;


%% Histogram of responses
%------------------------------------------
if enablehist
snridx = find(soisnrdb>=0,1);
erridx = find(abs(direrr_array)<1e-4,1); % No error
scenidx = 1;

figure;
v = histopt{snridx}.hist;
magbins = histopt{snridx}.magbins;
bar(magbins,v','histc');
[i,j] = find(v~=0);
j = max(j);
set(gca,'xlim',[magbins(1) magbins(j)])
title('Ideal MVDR')

figure;
v = histsphere{2,snridx,erridx}.hist;
magbins = histopt{snridx}.magbins;
bar(magbins,v','histc');
[i,j] = find(v~=0);
j = max(j);
set(gca,'xlim',[magbins(1) magbins(j)])
title('Multiple constraints')


figure;
v = histsphere{1,snridx,erridx}.hist;
magbins = histopt{snridx}.magbins;
bar(magbins,v','histc');
[i,j] = find(v~=0);
j = max(j);
set(gca,'xlim',[magbins(1) magbins(j)])
title('Outer constraint only')

% % Response in ring (essentially response over outer hyperball if n>>1)
% figure;
% vr = v(2,:) - (con{2}.radius(1)/con{2}.radius(2))^numants*v(1,:);
% bar(magbins,[vr' v(2,:)'],'histc');
% set(gca,'xlim',[magbins(1) magbins(j)])
% title('Outer constraint only - Ring response')

figure;
v = histsphere{3,snridx,erridx}.hist;
magbins = histopt{snridx}.magbins;
bar(magbins,v','histc');
[i,j] = find(v~=0);
j = max(j);
set(gca,'xlim',[magbins(1) magbins(j)])
title('Inner constraint only')

figure;
v = histmvdr{snridx,erridx}.hist;
magbins = histopt{snridx}.magbins;
bar(magbins,v','histc');
[i,j] = find(v~=0);
j = max(j);
set(gca,'xlim',[magbins(1) magbins(j)])
title('MVDR with steering error')
end