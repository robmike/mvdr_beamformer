function bcast_test()
%% These examples from http://www.scipy.org/EricsBroadcastingDoc
    bcast(0:2,[0:10:30]','+')
    bcast([0:10:30]',0:2,'+')
%%   
    codes = [225, 80;
            290,76;
            100,61;
            125,68];
    obs = [245,74];
        
    d = bcast(codes,obs,'-');
    dist = sqrt(sqrt(sum(d.^2,2)));
    [t,idx] = min(dist)
    % idx should be 1
    
%% Next example based on Figure 6 but, since no numbers are provided, we
% make them up ourselves.
    codes = [225, 80;
            290,76;
            100,61;
            125,68];
        
    obs = [ 245,74;
            140,65;
            140,75;
            260,75];
        
     % Augment dimensions of arrays
    [m,n]=size(codes);
    rcodes = zeros([m,1,n]);
    rcodes(:,1,:) = codes;
    %[m,n]=size(obs);
    %robs = zeros([1 m n]);
    %robs(1,:,:) = obs;
    robs = shiftdim(obs,-1);
    d = bcast(rcodes,robs,'-');
    dist = sqrt(sqrt(sum(d.^2,3)));
    [t,idx] = min(dist) % Should be 1 4 4 2
%%
end