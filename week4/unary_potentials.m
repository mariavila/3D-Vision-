function [ unary_pot ] = unary_potentials( left_image,right_image, min_disparity, max_disparity, window_size)
    [left_rows, left_cols] = size(left_image);
    unary_pot = zeros(size(left_image,1) * size(left_image,2), max_disparity+1);
    
    % we need to pad half the window size 
    padding = floor(window_size/2); 
    
    left_image = padarray(left_image,[padding padding]);
    right_image = padarray(right_image,[padding padding]);
    
    % init the weights
    weights = zeros(window_size);
    weights(:) = 1/(window_size * window_size);

    for row = 1+padding:left_rows+padding
        for col = 1+padding:left_cols+padding
            pixel_idx = (col-1)*left_rows+row;
            % get patch from with current pixel as center
            window_left = double(left_image(row-padding:row+padding,col-padding:col+padding));            
            
            % calculate the indexes for sliding window
            %min_col = max(max(1+padding, min_disparity), col - max_disparity);
            min_col = max(1+padding, col + min_disparity);
            max_col = min(left_cols + padding, col + max_disparity);
            
            % init the cost to the worst for each matching algorithm
            disparity = zeros(1,max_disparity+1);
            disparity_assigned = false(1,max_disparity+1);
            for kk = min_col:max_col
                window_right = double(right_image(row-padding:row+padding,kk-padding:kk+padding));
                
                % compute cost for every pixel of the window
                ssd = sum(sum(weights.*double((window_left-window_right).^2)));
                disparity(kk - min_col+1) = ssd;
                disparity_assigned(kk - min_col +1) = true;
            end
            %Assign disparity possibility          
            max_d = sum(disparity(disparity_assigned));
            unary_pot_aux = zeros(1,max_disparity+1);
            unary_pot_aux(disparity_assigned) = (max_d - disparity(disparity_assigned)) / max_d; 
            unary_pot(pixel_idx, :) = unary_pot_aux;
        end        
    end
end