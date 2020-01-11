function [F, idx_inliers] = ransac_fundamental_matrix(x1, x2, threshold)

    [Ncoords, Npoints] = size(x1);
    max_it = 1000;

    % ransac
    it = 0;
    best_inliers = [];
    % probability that at least one random sample set is free of outliers
    p = 0.999; 
    while it < max_it

        points = randomsample(Npoints, 8);
        F = fundamental_matrix(x1(:,points), x2(:,points));
        F = F / norm(F);
        inliers = compute_inliers(F, x1, x2, threshold);

        % test if it is the best model so far
        if length(inliers) > length(best_inliers)
            best_inliers = inliers;
        end    

        % update estimate of max_it (the number of trials) to ensure we pick, 
        % with probability p, an initial data set with no outliers
        fracinliers =  length(inliers)/Npoints;
        pNoOutliers = 1 -  fracinliers^8;
        pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
        pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
        max_it = log(1-p)/log(pNoOutliers);

        it = it + 1;
    end

    % compute H from all the inliers
    F = fundamental_matrix(x1(:,best_inliers), x2(:,best_inliers));
    F = F / norm(F);
    idx_inliers = best_inliers;

end

function idx_inliers = compute_inliers(F, x1, x2, threshold)
    x1 = normalise(x1);
    x2 = normalise(x2);
    
    l1 = F' * x2;
    l2 = F * x1;
    
    l1 = l1 ./ sqrt(repmat(l1(1,:).^2 + l1(2,:).^2 + l2(1,:).^2 + l2(2,:).^2, size(l1,1), 1));
    l2 = l2 ./ sqrt(repmat(l1(1,:).^2 + l1(2,:).^2 + l2(1,:).^2 + l2(2,:).^2, size(l1,1), 1));
    
    d1 = abs(sum(l1.*x1));
    d2 = abs(sum(l2.*x2));
    d = d1 + d2; 
    
    idx_inliers = find(d < threshold);

end

function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

end
    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat
end