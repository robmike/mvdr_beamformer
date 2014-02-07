function z = normfunc(theta,phi,numants)
    
    M = numants;

    if ~isvector(theta) || ~isvector(phi)
        error('Inputs must be vectors')
    end

    if (isrow(theta) && isrow(phi)) || (iscol(theta) && iscol(phi))
        phi = phi';
    end
    % At this point one vector is row and one vector is column

    ed = pi*bcast(sin(theta),sin(phi),'-');
    %z = cos((M-1)*ed/2) .* sin( M * ed /2 )./sin(ed/2);
    %z = sqrt(2*(numants-z));
    
    % Alternatively
    z = sqrt(2*M-1-(2*M-1)*diric(ed,(2*M-1)));