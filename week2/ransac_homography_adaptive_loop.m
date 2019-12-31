function [H, idx_inliers] = ransac_homography_adaptive_loop(x1, x2, th, max_it)

[Ncoords, Npoints] = size(x1);

% ransac
it = 0;
best_inliers = [];
while it < max_it
    
    points = randomsample(Npoints, 4);
    H = homography2d_2nd(x1(:,points), x2(:,points)); % ToDo: you have to create this function
    inliers = compute_inliers(H, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    p=0.99;
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute H from all the inliers
H = homography2d(x1(:,best_inliers), x2(:,best_inliers));
idx_inliers = best_inliers;


function idx_inliers = compute_inliers(H, x1, x2, th)
    % Check that H is invertible
    if abs(log(cond(H))) > 15
        idx_inliers = [];
        return
    end
    

    % compute the symmetric geometric error
    n_points = size(x1);
    n_points = n_points(2);
    for i = 1:n_points
        x1p = H*x1(:,i);
        x2p = H*x2(:,i);
        d2_1 = (x1(1,i)/x1(3,i) - x1p(1)/x1p(3))^2+(x1(2,i)/x1(3,i) - x1p(2)/x1p(3))^2;
        d2_2 = (x2(1,i)/x2(3,i) - x2p(1)/x2p(3))^2+(x2(2,i)/x2(3,i) - x2p(2)/x2p(3))^2;
        d2(i) =(d2_1 + d2_2);  
    end
    idx_inliers = find(d2 < th.^2);


function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat

function homo = homography2d(x1,x2)
    n_points = size(x1);
    n_points = n_points(2);
    x1p = normalise(x1);
    x2p = normalise(x2);
    mean1 = mean(x1,2);
    mean2 = mean(x2,2);
    std1 = std(x1,1,2);
    std2 = std(x2,1,2);
    s1 = std1./sqrt(2);
    s2 = std2./sqrt(2);
    T1 = [s1(1), 0, -mean1(1);...
          0, s1(2), -mean1(2);...
          0, 0, 1];
    T2 = [s2(1), 0, -mean2(1);...
          0, s2(2), -mean2(2);...
          0, 0, 1];
    x1p = T1 * x1p;
    x2p = T2 * x2p;
    count = 1;
    for i = 1:n_points
        x1i = x1p(:,i);
        x2i = x2p(:,i);
        A(count,:) = [0, 0, 0,...
                     -x2i(3)*x1i(1), -x2i(3)*x1i(2), -x2i(3)*x1i(3), ...
                      x2i(2).*x1i(1), x2i(2).*x1i(2), x2i(2).*x1i(3)];
                  
        A(count+1,:) = [x2i(3)*x1i(1), x2i(3)*x1i(2), x2i(3)*x1i(3), ...
                        0, 0, 0,...
                        x2i(1)*x1i(1),x2i(1)*x1i(2),x2i(1)*x1i(3)];
        count = count + 2;
    end
    [U,D,V] = svd(A);
    h = V(9,:);
    homo = reshape(h,3,3);
    homo = inv(T2)*homo*T1;
    
function homo = homography2d_2nd(x1,x2)
    n_points = size(x1);
    n_points = n_points(2);
    %x1n = normalise(x1);
    %x2n = normalise(x2);
    x1n = x1;
    x2n = x2;
    count = 1;
    A = zeros(8,9);
    for i = 1:n_points
        x1i = x1n(:,i);
        x2i = x2n(:,i);
        A(count,:) = [0, 0, 0, -x2i(3).*x1i', x2i(2).*x1i'];
        A(count+1,:) = [x2i(3).*x1i', 0, 0, 0, -x2i(1).*x1i'];
        count = count + 2;
    end
    [~,~,V] = svd(A);
    h = V(:,9);
    homo = reshape(h,3,3);
    homo = homo';
    %homo = T2\homo*T1;