function [ diff ] = gs_errfunction( P0, Xobs )

    % First get the H matrix
    H = reshape(P0(1:9), [3,3]);
    % get x and xp contained in Xobs
    n_points = size(Xobs,1) / 2;
    x = Xobs(1:n_points); % x
    x = reshape(x, [2,size(x,1)/2]); % from nx1 to 2x(n/2)
    x = [x; ones(1, length(x))]; %Add third coordinate equal to 1
    xp = Xobs(n_points+1:end); % x'
    xp = reshape(xp, [2,size(xp,1)/2]); % from nx1 to 2x(n/2)
    xp = [xp; ones(1, length(xp))]; %Add third coordinate equal to 1
    % reprojection error
    xhat = P0(10:end);
    xhat = reshape(xhat, [2,size(xhat,1)/2]); % from nx1 to 2x(n/2)
    xhat = [xhat ; ones(1,size(xhat,2))]; % from euclidean to homogeneous
    xhatp = H*xhat;
    xhatp = xhatp./xhatp(3,:);
%     xhat = inv(H)*xhatp;
%     xhat = xhat./xhat(3,:);
    diff = sqrt(sum((x(1:2,:)-xhat(1:2,:)).^2,1) + sum((xp(1:2,:)-xhatp(1:2,:)).^2,1));
%     disp(sum(diff.^2))
end