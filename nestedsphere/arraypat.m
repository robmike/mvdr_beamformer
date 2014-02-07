function arraypat(w,varargin)
% Display the array pattern for beamformer with coefficients w
    res=0.01;
    numants = length(w);
    t = -pi:res:pi;
    
    x = bsteerdir(t',numants);;
    
    semilogy(t./pi,abs(w'*x),t./pi,ones(1,length(t)),':');
    xlabel('Rads/\pi')
    
    %% Plot vertical lines 
    if nargin>1
        lx = varargin{1}./pi;
        lx = torow(lx);
        n = length(lx);
        y = get(gca,'ylim')';
        line(repmat(lx,2,1),repmat(y,1,n),'linestyle',':');
    end
    
end