function y = bsteer(theta,numants)

    m = repmat((0:numants-1)',1,length(theta));
    theta = repmat(theta,numants,1);
    %theta = pi*sin(theta);
    y = exp(j*m.*theta);