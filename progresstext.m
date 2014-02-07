function [timeleft ratios] = progresstext(loopidx,loopsize,varargin)
% [timeleft ratios] = progresstext(loopidx,loopsize,varargin)
% Print the progress and estimated time remaining for each loop
% loopidx
% Entries must be ordered from innermost loop to outermost
% Example: progresstext([innerloopidx outerloopidx],[inloopsize
% outloopsize])


    persistent starttimes last_t lastratios;
    teps = 1; % Limit no. of updates to every teps seconds
    fid = 1; % Screen

    yesprint = true;
    if nargin>2
        yesprint = varargin{1};
    end

    n = length(loopidx);

    % Calculate percentage completion of each loop
    % Computes linear index of each set of n innermost loops divided by size
    % of product of these loops, for n=1:numloops
    ratios = cumsum([loopidx(1) (loopidx(2:end)-1).*cumprod(loopsize(1:end-1))]);
    ratios = ratios./cumprod(loopsize);


    if isempty(starttimes)
        starttimes = repmat(clock,n,1);
    end
    if isempty(last_t)
        last_t = clock; % Will cause early return for very first iteration
    end
    if isempty(lastratios)
        lastratios = zeros(1,n);
    end
    
  
    % Reset/initialize
    if all(ratios==0) || n~=length(lastratios)
        starttimes = repmat(clock,n,1);
        lastratios = zeros(1,n);
        last_t = clock;
% $$$     elseif any(ratios<lastratios) % Counter wrap around
% $$$         zidx = find(ratios<lastratios); % Counter wrap around
% $$$         starttimes(ratios<lastratios,:)=repmat(clock,length(zidx),1);
        end



    if etime(clock,last_t) < teps % Do not update too frequently
        return;
    end


    ratios(ratios>1) = 1;
    ratios(ratios<0) = 0;
    
%if all(ratios>=1)
%    return;
%end


% FIXME: When an inner loop wraps around we should use time
% difference from last update and lastratio to calculate estimated
% time remaining so that it does not give an inaccurate value on
% the first iteration as it does now.

 p = ratios - lastratios;
p(p<0) = 1 + p(p<0);
p(p<1e-5)=0.01; % Avoid division by zero

timeleft = etime(clock,last_t)*(1-ratios)./p;
%     timeleft=zeros(1,n);
%nzidx = find(ratios~=0);

% $$$     for i = nzidx
% $$$         runtime = etime(clock,starttimes(i,:));
% $$$         timeleft(i) = runtime/ratios(i) - runtime; 
% $$$ 
% $$$         % Todo: Use a FIR filter to better estimate remaining time
% $$$     end

    lastratios = ratios;
    last_t = clock;

    % Print time left, outermost loop prints leftmost
    if yesprint
        printtimeleft(fliplr(timeleft),fliplr(ratios),fid);
    end


end

function printtimeleft(timeleft,ratios,fid)
% Prints percentage complete and estimated time remaining
    for i = 1:length(timeleft)
        fprintf(fid,'%s (%.0f%%) \t',sec2timestr(timeleft(i)),100*ratios(i));
    end
   fprintf(fid,'\n');
end


function tstr = sec2timestr(sec)
% Adapted from progressbar.m by Steve Hoelzer (no longer contains
% much of that code)
% Convert a time measurement from seconds into a human readable string.

% Converts number of seconds into days,hours,minutes,seconds
d = floor(sec/86400); % Days
sec = sec - d*86400;
h = floor(sec/3600); % Hours
sec = sec - h*3600;
m = floor(sec/60); % Minutes
sec = sec - m*60;
s = floor(sec); % Seconds


% Maybe this can made more Matlabish by removing the if statements by
% doing something like this:
% v = [d h m s]
% idx = findfirst(v); % You would need to write a fast function to
% find first non-zero entry (This is the main problem with this method)
% tstr = [sprintf('%i:',v(idx:end-1) sprintf('%i',s)]

if d>0
    tstr = sprintf('%02id %02ihr',d,h);
elseif h>0
    tstr = sprintf('%02i:%02i:%02i',h,m,s);
else
    tstr = sprintf('%02i:%02i',m,s);
end

end
