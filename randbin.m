function r=randbin(m,n)
% Streamlined implementation of specific case of randint
% Generates m by n matrix of equiprobable elements in the set {0,1}
r = floor(2*rand(m, n)); % Floor is faster than round I believe