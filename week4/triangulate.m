function [ X ] = triangulate( x1, x2, P1, P2, imsize )
% Performs a triangulation with the homogeneous algebraic method (DLT)
%    The entries are (x1, x2, P1, P2, imsize), where:
%        - x1, and x2 are the Euclidean coordinates of two matching 
%           points in two different images.
%        - P1 and P2 are the two camera matrices
%        - imsize is a two-dimensional vector with the image size
    nx = imsize(1);
    ny = imsize(2);
    H = [2/nx 0 -1; 0 2/ny -1; 0 0 1];
    x1 = [x1(1); x1(2); 1];
    x2 = [x2(1); x2(2); 1];
    x1s = H * x1;
    x2s = H * x2;
    P1s = H * P1;
    P2s = H * P2;
    A = [x1s(1) * P1s(3, :) - P1s(1, :);
         x1s(2) * P1s(3, :) - P1s(2, :);
         x2s(1) * P2s(3, :) - P2s(1, :);
         x2s(2) * P2s(3, :) - P2s(2, :)];
    [~,~,V] = svd(A);
    X = V(:,4);
    X = X./X(end);
end

