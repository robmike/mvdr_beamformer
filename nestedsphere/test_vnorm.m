% Tests vnorm.m for 2-dimensional case only
function test_vnorm()
n=10;
A = 10*rand(3,10)-5;

sizA = size(A);
dims = 1:length(size(A));
type = {inf,-inf,1,2,3,4};

disp('Testing...');

for d=dims
   for t = [type{:}]
      issame(A,d,t);
   end
end

A =  [8 1 6; 3 5 7; 4 -9 2];

%%
for i = 1:size(A,2)
   norm(A(:,i),2)
end
%%


disp('No errors')

end

% -----------------------------------------
function issame(A,dim,type,varargin)

if(dim==2)
    A=A';
end

if ~isempty(varargin)
   x = varargin{:}; 
else
x = vnorm(A,1,type);
end

y = zeros(1,size(A,2));
for i = 1:size(A,2)
    y(i) = norm(A(:,i),type);
end

if any( abs(x-y)>1e-9)
   error('Bug in code')
   fprintf('%f ',x-y)
   fprintf('\n');
end


end
