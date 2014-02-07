function coeffs = socp_coeffs(presumed,atten,corrmat,ellipses,varargin)
% coeffs = socp_coeffs(presumed,atten,corrmat,ellipsemat)
% 'presumed' is matrix with ith column being the center (presumed
% steering vector) of the  ith ellipse.
% ellipsemat is 3-dimensional array. (:,:,i) holds the matrix
% describing the ith ellipse

logfile=0;
%logfile=fopen('socp.log','w');
%if(nargin>4)
%    logfile = varargin{1};
%end

% If imagzero is true the imaginary part of the inner product
% . This results from a restriction on non-attenuation of the
% magnitude of the response in the ellipse. Without this constraint
% only the real part of of the response is restricted to be 
% 'non-attenuated'.
%imagzero = true;

n = size(corrmat,1);


if(size(presumed,1)==1 && size(presumed,2)~=1)
    error('Presumed vector must be a column vector');
end


% If ellipses is a vector, assume that length(ellipses) spheres are 
% desired with radius given by entries in vector
if(isvector(ellipses))
    numellipse = length(ellipses);
    ellipsemat = zeros([n n numellipse]); 
    for i = 1:numellipse
        ellipsemat(:,:,i) = ellipses(i)*eye(n); 
    end
    
    if( norm(presumed) - max(ellipses) < 0.1 )
        fprintf(['Hyperplane does not appear to ' ...
                            'intersect all of the cones. ' ...
                            'Problem seems infeasible. Continuing ' ...
                            'anyway\n'])
        fprintf('Max cone grad: %f, Hyperplane grad:%f\n\n', ...
                max(ellipses),norm(presumed));
    end

else
    ellipsemat = ellipses;
    numellipse = size(ellipses,3);
end



if(length(atten)~=numellipse)
    error('Expected %i matrices describing ellipses, received %i'...
              ,length(atten),numellipse);
end

numrows = (n+1)*(numellipse+1); % n+1 for each cone constraint

At = zeros(numrows,n+1);
c = zeros(numrows,1); 
K.q = n+1; K.q = K.q(ones(1,numellipse+1)); % Repeat n+1
                                           % numellipse+1 times

%    K.xcomplex = 1:numrows; 
%    K.ycomplex = 2:n; 


U = chol(corrmat);

c(n+2:n+1:end) = -atten;
At(1,1)=-1;
At(2:n+1,2:n+1) = -U;

%size(At(n+2:n+1:end,2:end))
%size(repmat(-presumed',numellipse,1)) 
At(n+2:n+1:end,2:end) = -repmat(presumed,1,numellipse)';


ridx = n+3;
for i = 1:numellipse
    At(ridx:ridx+n-1,2:end) = -ellipsemat(:,:,i).';
    ridx = ridx + n+1;
end

b = [-1 ; zeros(n,1)]; % Maximize first element


if(exist('imagzero','var') && imagzero)
    c=[0;c];
    At = [0 -[presumed(n/2+1:end)' -presumed(1:n/2)'] ; At];
    K.f = 1; % Restricts first row of At to be zero
    % K.xcomplex = 1:(numrows+1)
end

%%% Solve the second order cone programming problem using Sedumi
%% min tau s.t. norm(Ux)<=tau,
%norm(A1'x)<=c'x-atten(1),norm(A2'x)<=c'x-atten(2)
%% c=presumed, A(i) = varargin{i}

% Sedumi options
pars=[];
pars.fid=logfile; % 0 for quiet mode (no output to screen),1 for screen
% Algorithms for pars.alg: 0=1st order wide region algorithm, 2 =
% xz-predictor corrector, 1=centering predictor-corrector
%pars.alg=0; % 0 is the only algorithm that seems to work
pars.eps=1e-8; % Accuracy


% In Sedumi problem is set in the form:
% max b'*y such that c - A'y is in K, where K is a (set of) second
% order cones
[x,y,optinfo] = sedumi(At,b,c,K,pars);
% Sedumi automatically detects correct orientation of A (as long as
% it is not square) so we can pass it either A or A'. 
if(any(eigK(c-At*y,K) < -0.1)) %|| any(eigK(x,K) < -0.1))
    warningtext = 'eigK feasibility check failed\n';
    fprintf(warningtext);
    %warning('SOCP:feascheck',warningtext);
    fprintf('%f ',eigK(c-At*y,K)'); fprintf('\n');
    %fprintf('%f ',eigK(x,K)'); fprintf('\n');
end

%fprintf('Feasratio=%f\n',optinfo.feasratio);
if(optinfo.feasratio<0.8)
    %warningtext = sprintf('Feasibility ratio is %f\n',optinfo.feasratio);
    %warning('SOCP:feasratio',warningtext);
    %fprintf(warningtext);
end

dualgap = c'*x - b'*y;
if(dualgap > 0.1 || dualgap<-0.01) % Should be positive, ideally should be zero
    warningtext = sprintf('Dual gap is %f\n',dualgap);
    %warning('SOCP:dualgap',warningtext);
    fprintf(warningtext);
end


% If dual is infeasible or other convergence error
if(optinfo.numerr==2 || optinfo.dinf==1) 
    warning('SOCP:noconverge', ['Coefficient optimization did not ' ...
                        'converge'])
end

%% Another implementation using fmincon (doesn't need Sedumi, not
%necessary to use them both simultaneously aside from testing/debugging).
confunc = @(z) socpcon(z,c,At,K);

options = [];
options = optimset(options,'TolFun',1e-4,'MaxFunEvals',400,'MaxIter',400);
%options = optimset(options,'Display','notify');
%x0 = [0;presumed(:,1)];
%x0(1) = norm(x0);
x0 = zeros(n+1,1);


%[y2 v exitflag optinfo2]=fmincon(@(x) [1;zeros(n,1)]'*x, x0,[],[],[],[],[],[],confunc,options);
%if(norm(y2-y)>0.1)
%    warning('fmincon and sedumi do not give same results');
%end


if(logfile~=0 && logfile~=1)
    fclose(logfile);
end


coeffs = y(2:end);
