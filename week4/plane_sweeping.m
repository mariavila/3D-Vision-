function [ disparity ] = plane_sweeping(I1, I2, P1, P2, range_depth, size_window, cost_function, step_depth)


r_window = ceil(size_window/2) - 1;
sampling_depth = range_depth(1):step_depth:range_depth(2);
[heigth, width] = size(I2);

% In third dimension, first value is for disparity/depth and second for min
% cost
switch cost_function
    case 'SSD'
        disparity_computation = Inf*ones([heigth width 2]);
    case 'NCC'
        disparity_computation = -Inf*ones([heigth width 2]);   
end
for k = 1:length(sampling_depth)
    % Compute the homography that relates the two images
    kk = sampling_depth(k);
    PI = P1(3, :) - [0 0 0 kk];
    A = pinv([P1; PI]);
    A_hat = A(:, 1:3);
    H = P2*A_hat;
    % Take into account those points that fall inside the area of the other
    % image
    corners = [1 width 1 heigth]; 
    I_reprojected = apply_H_v2(I2, pinv(H), corners);

    %For each point in the image
    for i = 1:heigth
        for j = 1:width
            window_left = I1(max(i - r_window, 1 + r_window):...
                min(i + r_window, heigth - r_window),...
            max(j - r_window, 1 + r_window):...
            min(j + r_window, width - r_window));
        
            window_right  = I_reprojected(max(i - r_window, 1 + r_window):...
                min(i + r_window, heigth - r_window),...
            max(j - r_window, 1 + r_window):...
            min(j + r_window, width - r_window));
        
            size_block = size(window_right);
            
            weights = (1/(sum(size_block)))*ones(size_block);
            switch cost_function
                case 'SSD'
                    ssd = sum(sum(weights.*double((window_left-window_right).^2)));
                        if ssd < disparity_computation(i, j, 2)
                            disparity_computation(i, j, 1) = kk;
                            disparity_computation(i, j, 2) = ssd;
                        end
                case 'NCC'
                    sum_left = sum(window_left(:).*weights(:));
                    sum_right = sum(window_right(:).*weights(:));

                    sigma_left = sqrt(sum( weights(:).* (window_left(:) - sum_left).^2 ));
                    sigma_right = sqrt(sum( weights(:).* (window_right(:) - sum_right).^2 ));

                    ncc = sum( weights(:).*(window_left(:)-sum_left).*(window_right(:)-sum_right) )/(sigma_left*sigma_right);
                   
                    if ncc > disparity_computation(i, j, 2)
                        disparity_computation(i, j, 1) = kk;
                        disparity_computation(i, j, 2) = ncc;
                    end
            end            
        end
    end
end

disparity = disparity_computation(:, :, 1);

end