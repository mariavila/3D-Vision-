function [ error ] = gs_errfunction(P, x, X)
    
    x_hat = euclid(P * X);
       
    error = sqrt(sum((x - x_hat) .^ 2));
end