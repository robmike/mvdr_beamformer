function [val,idx]=findclosest(x,target,varargin)
    % Find value in array x that is closest (in 1-norm) to scalar 'target'.
    % Returns (linear) index of closest match. 
    % Optional argument specifying the largest mismatch allowed
    
    e = inf;
    if nargin>2
       e = varargin{1}; % maximum 
    end
    
    if ~isscalar(target)
        error('Target must be scalar')
    end

    [d,idx] = min(abs(x(:)-target));
    if d > e % If match is not within constraints
        val=[];
        idx=[];
    else
        val = x(idx);
    end
    
end