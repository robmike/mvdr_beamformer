m = 1e6;
n=2;

rm = .5;

r = sqrt(rand(1,m))*rm;
u = 2*pi*rand(1,m);
[x,y]=pol2cart(u,r);

ve=-rm:.1:rm;

mhist2d = hist2d([x;y]',ve,ve);
plot2dhist(mhist2d,ve,ve,'x','y','example');  colorbar

return
%%

    mYX = rand(10000,2);
    vXEdge = linspace(0,1,100);
    vYEdge = linspace(0,1,200);
    mHist2d = hist2d(mYX,vYEdge,vXEdge);
 
    Plot2dHist(mHist2d, vXEdge, vYEdge, 'X', 'Y', 'Example'); colorbar