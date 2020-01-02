function [ diff ] = gs_errfunction( P0, Xobs )

    % First get the H matrix
    H = reshape(P0(1:9), [3,3]);

    % get x1 and x2 contained in Xobs
    n_points = size(Xobs,1) / 2;
    x1 = Xobs(1:n_points); % x
    x1 = reshape(x1, [2,size(x1,1)/2]); % from nx1 to 2x(n/2)
    x2 = Xobs(n_points+1:end); % x'
    x2 = reshape(x2, [2,size(x2,1)/2]); % from nx1 to 2x(n/2)
    
    % reprojection error
    xhat = P0(9+1:end);
    xhat = reshape(xhat, [2,size(xhat,1)/2]); % from nx1 to 2x(n/2)
    xhat = [xhat ; ones(1,size(xhat,2))]; % from euclidean to homogeneous
    xhatp = H*xhat;
    diff = sum(x1-euclid(xhat)).^2 + sum(x2-euclid(xhatp)).^2;

end