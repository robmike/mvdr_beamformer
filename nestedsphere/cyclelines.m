function cyclelines(h)
% CYCLELINES - Adjust the line properties of a plot
%   h is a figure handle
% FIXME: Construct all unique permutations of dashes and markers
% FIXME: Marker size changes thickness when changing line thickness
dash = {'--'}; %{'--','-','-.',':'};
lthick = 2; % FIXME: Get current default
colour = 'brcmk'; %rgbcmyk
marker = '+oxds^><*'; %{'+ox*sd^v><ph'}
msize = 8; % FIXME: Get current default

cprod = setprod(1:length(dash),1:length(marker));
ncombos = length(cprod);

lns = findobj(h,'type','line');
legh = findobj(findobj(h,'tag','legend'),'type','line');
lns = setdiff(lns,legh);

for i=1:length(lns)
    %didx = rem(i,length(dash))+1;
    %midx = rem(i,length(marker))+1;
    cidx= rem(i-1,length(colour))+1;
    modidx = rem(i-1,ncombos)+1;
    didx = cprod(modidx,1);
    midx = cprod(modidx,2);

    set(lns(i),'LineStyle',dash{didx},'marker',marker(midx), ...
               'color',colour(cidx),'markersize',msize,'linewidth',lthick);
    %set(lns(i),'Marker','None');
end
