function [ diff ] = gs_errfunction( P0, Xobs )

    % First get the H matrix
    H = reshape(P0(1:9), [3,3]);
    % get x1 and x2 contained in Xobs
    n_points = size(Xobs,1) / 2;
    x1 = Xobs(1:n_points); % x
    x1 = reshape(x1, [2,size(x1,1)/2]); % from nx1 to 2x(n/2)
    x1 = [x1; ones(1, length(x1))];
    x2 = Xobs(n_points+1:end); % x'
    x2 = reshape(x2, [2,size(x2,1)/2]); % from nx1 to 2x(n/2)
    x2 = [x2; ones(1, length(x2))];
    % reprojection error
    xhat = H*x1;
    xhat = xhat./xhat(3,:);
    xhatp = inv(H)*x2;
    xhatp = xhatp./xhatp(3,:);
    diff = sqrt(sum((x1(1:2,:)-xhatp(1:2,:)).^2,1) + sum((x2(1:2,:)-xhat(1:2,:)).^2,1));
end