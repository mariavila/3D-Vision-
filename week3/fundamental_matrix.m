function [F_es] = fundamental_matrix(x1, x2)
    [x1Norm, t1] = normalise2dpts(x1); 
    [x2Norm, t2] = normalise2dpts(x2); 
    x1 = x1Norm(1,:)';
    y1 = x1Norm(2,:)';
   
    x2 = x2Norm(1,:)';
    y2 = x2Norm(2,:)';

    W = [x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(size(x1,1),1)];
    
    [~,~, V] = svd(W);

    V = V(:,9);
    F_es = reshape(V, 3,3)';
    [U, D, V] = svd(F_es);
    D(3,3) = 0;
    F_es = U * D * V';
            
    F_es = t2' * F_es * t1;
end
