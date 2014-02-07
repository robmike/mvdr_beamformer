%function polishfigs()
% Make figs publication ready.
epsexport = true;
plotdir ='plots/';
figwidth = 8.89; % 3.5 inches 

xfontsize = 14;
yfontsize = xfontsize;
tickfontsize = 14;

% Make legend size really big if using LaPrint so that the Legend box is
% big enough
lfontsize = 14; %Legend text size


dash = {'--','-','-.',':'};
lthick = 2; % FIXME: Get current default
colour = 'k'; %rgbcmyk

% Create plot direcotory if it does not exist
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end

fprintf('Polishing figs...\n');

figs=findobj(0,'type','figure'); % Get all figures
    
for f=torow(figs)
        axl = findobj(f,'type','axes');
    leg = findobj(axl,'tag','legend');
    ax = setdiff(axl,leg); % Axes which are not legends
 
       lns = findobj(f,'type','line');
   
% ---- Remove markers
            set(lns,'Marker','None');
    

% ---- Convert lines to dash-dot and black and white
    set(lns,'color','k','linewidth',lthick);
    
i=0;
for l=torow(lns)
    % ---- We have only 4 line types in Matlab. Convert ideal MVDR to pseudo
% dashed line by using blank lines and markers.
    try
    dname = get(l,'displayname');
    catch
        dname = '';
    end
    switch dname
        case legtext{end} % Is it the ideal MVDR?
            set(l,'linestyle','none','Marker','.','markersize',8); 
        case '2-constraint' % MVDR with mismatch
            set(l,'linestyle','-.','Marker','+','markersize',8);
        otherwise
            set(l,'linestyle',dash{mod(i,length(dash))+1});
            i=i+1;
    end
end
    
    
    for a=torow(ax)
    % ---- Increase font size on axis
          h = get(a,'xlabel');
          set(h,'fontsize',xfontsize);
       
          h = get(a,'ylabel');
          set(h,'fontsize',yfontsize);
          
          set(a,'fontsize',tickfontsize); % Why does this change legend font size?
          
          h = legend(a);
    set(h,'fontsize',lfontsize);
    end
    
%     legtext = findobj(get(leg,'children'),'type','text');
%     olfsize = max(cell2mat(get(legtext,'fontsize')));
%     
%     if olfsize~=0
%     ratio = lfontsize/olfsize;
%     %set(legtext,'fontsize',lfontsize);
%     % We need to resize legend box at the same time
%     resize_legend(leg,ratio);
%     end
    
    
     % ---- Remove titles (actually just hide them)
     t = get(a,'title');
     set(t,'visible','off')
    
 




    % (Courtesy of P. Kabal Matlab plotting routines)
    % To avoid lines in the legend being redrawn when printing
    %  - Turn off the legend propery tag
    %  - Make the legend axes handle invisible
    % Leave the legend on top so that it is not hidden by the plot
    %set (leg, 'Tag', '');
    %set (leg, 'HandleVisibility', 'off');



%----- Export to eps
if epsexport
   filename = get(f,'userdata'); 
   print(f,'-depsc2',[plotdir filename '.eps']) 
end


%----- Print figure using laprint
% laprint(f,filename,...
%     'width',8.89,... % width of image in latex document
%     'factor',0.8,... % Ratio between latex doc image and eps
% %    'scalefonts','on',...
%     'keepfontprops','on',... % Use Matlab font properties?
%     'createview','on',... % Automatically produce a latex document
%     'processview','on',... % Compile said latex document
%     );
end


fprintf('Done polishing.\n');
%end % End function
