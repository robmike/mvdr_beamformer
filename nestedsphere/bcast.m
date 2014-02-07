function z = bcast(x,y,varargin)
% BCAST - Rudimentary implementation of SciPy's 'broadcasting': Binary
% operation applied to 'intuitive' repmat-ting of arguments.
%   BCAST(X,Y,OP) applies the binary operation OP to repmat-ted X and Y
%   arguments. The repmat-ting is done so as to match the dimensions of the
%   transformed X and Y in the obvious way. (Explained futher below).
%
%   The corresponding dimensions of X and Y must be equal or else one of
%   them singleton (The smaller dimension vector has leading ones padded
%   to make both 'size' vectors the same length). See below for examples.
%
%   OP may be a function handle, or one of the following strings (separated
%   by commas):
%   +,-,*,.*,^,.^,\,/,./,.\
%   (Type 'help +' or 'help plus' for more info about these operations)
%
%   For the first calling signature, OP must have a signature OP(X,Y) or 
%   OP(X,Y,varargin)
%
%   BCAST(X,Y) and BCAST(X,Y,[]) are the same as BCAST(X,Y,'+')
%
%   BCAST(X,Y,OP,ARGOPS) does the same as above but passes the options 
%   ARGOPS to OP along with transformed X and Y
%   In this case OP must have a signature OP(X,Y,varargin) or
%   OP(X,Y,A1,A2,...,AN) where A1,...,AN are N mandatory arguments.
%
%   How is the tiling done?
%   The arguments are tiled (i.e. repmat-ted) such that both
%   array are scaled in their singleton dimensions as many times as the
%   corresponding dimension in the *other* argument. Thus, the resulting
%   transformed arguments have the same dimension (We
%   constrain non-singleton dimensions to be equal a priori).
%
%   WARNING: Matlab does not permit array with more than 2 dimensions to
%   have a trailing dimension of 1. Thus you may not pass, for example, 
%   an [M N 1] array to bcast. To work around this you can extend the
%   trailing dimension to the appropriate size yourself or permute the
%   dimensions of both arguments to move any trailing singleton dimensions
%   to a non-trailing location. See example below.
% 
%   EXAMPLE:    
%       % Matlab does broadcasting automatically for scalars
%       bcast(magic(4),3,'+') = magic(4) + repmat(3,4,4) = magic(4) + 3;
%  
%       % size(x) = [1 2]; size(y) = [3 1];
%       % size(z) = max([x;y],[],1) = [3 2];
%       z = bcast([1 2],[5:7]') % = [6 7; 7 8; 8 9]
%       
%       % size(x) = [2 2]; size(y) = [3 1];
%       % Yields an error since first dimension is not the same and neither
%       % is one. size(x) = [1 2 2] would also yield an error
%
%       % size(x) = [1 2] => [1 1 2]; % Ones padded dimension
%       % size(y) = [3 3 1];
%       % The desired outcome is that both x and y are repmat-ted along their singleton 
%       % dimensions resulting in [3 3 2] sized arrays. However, we cannot
%       % construct y with size [3 3 1] in Matlab. The solution is to pre-tile y:
%       y = repmat(y,[]1,1,2);
%       % We can't permute x and y in this case because then x would have a
%       % trailing singleton dimension. If size(x) was [3 2]
%       % then we could do permute the dimensions to make the singleton dimension
%       % non-trailing.
%       y = shiftdim(y,-1); 
%       x = shiftdim(shiftdim(x,-1),2);
%       z = bcast(x,y)
%       % size(z) = [2 3 3];
%       z = shiftdim(z,1); % Undo permutation
%
%       % Example of function handle and options
%       bcast(1:2,[1:3]',@(x,y,c) x+c.*y,3) % = [4 5; 7 8; 10 11]
%
%   See http://www.scipy.org/EricsBroadcastingDoc for better examples and
%   explanation including helpful diagrams.
%
%   This is a convenience function. It's speed may be suboptimal and may
%   require large amounts of memory. 
%
%   See also REPMAT PLUS


defaultop = @plus;

sx = size(x);
sy = size(y);

builtinop = false;

if nargin<=2
    op = defaultop;
    builtinop = true;
else
    op = varargin{1};
    varargin(1) = [];
    if isempty(op)
        op = defaultop;
    elseif ~isa(op,'function_handle')
        builtinop=true;
         % FIXME: These pairs should go in a dictionary (or whatever
         % Matlab's equivalent is).
        switch op
            case '+'
                op = @plus;
            case '-'
                op = @minus;
            case '*'
                op = @mtimes;
            case '.*'
                op = @times;
            case '^'
                op = @mpower;
            case '.^'
                op = @power;
            case '\'
                op = @mldivide;
            case '/'
                op = @mrdivide;
            case '.\'
                op = @ldivide;
            case './'
                op = @rdivide;
            otherwise
                error(['Operator must be a function handle or string '...
                    'representing a built-in arithmetic operator(help plus)']);
        end
    end
end

% Make sure (to best of our ability) that number of arguments is compatible
nops = 2+length(varargin); % Number of arguments to 
if nargin(op)>=0 && nargin(op)~=nops
    error('%s takes %i arguments. Tried to pass %i',func2str(op),nargin(op),nops);
end

% Built-in operators can broadcast themselves for scalars.
% All operators are can handle two scalars -- no broadcasting needed.
% FIXME: Is this an okay optimization for *all* operators?
if ( builtinop && (isscalar(x) || isscalar(y)) ) || (isscalar(x) && isscalar(y))
    z = op(x,y,varargin{:});
    return 
end

lsx = length(sx);
lsy = length(sy);


% FIXME: We can probably get rid of this 'if' and just do the padding for
% both. A negative padding will result in a no-op.
d = lsy - lsx;
if lsx<lsy
sx = [ones(1,d) sx]; % pad to make lsx == lsy
x = shiftdim(x,-d);
else
sy = [ones(1,-d) sy];
y = shiftdim(y,d);
end

% The docs say only the trailing dimension has to be equal or one but
% their examples imply all of them must.
if any(sx~=sy & sx~=1 & sy~=1)
    error('Corresponding dimensions must be the same, or one of them 1')
end



sz = max([sx;sy],[],1); % Dimension of result

z = zeros(sz); % Pre-allocate space for result

rx = repmat(x,sz-sx+1);
ry = repmat(y,sz-sy+1);


z = op(rx,ry,varargin{:});

