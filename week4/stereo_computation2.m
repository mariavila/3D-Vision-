function [disparity] = stereo_computation2(right_image,left_image, min_disparity, max_disparity, window_size, matching_cost)
    % Write a function called 'stereo_computation' that computes the disparity
    % between a pair of rectified images using a local method based on a matching cost 
    % between two local windows.
    % 
    % The input parameters are 5:
    % - left image
    % - right image
    % - minimum disparity
    % - maximum disparity
    % - window size (e.g. a value of 3 indicates a 3x3 window)
    % - matching cost (the user may able to choose between SSD and NCC costs)
    disparity = zeros(size(right_image));
    [left_rows, left_cols] = size(right_image);
    
    % we need to pad half the window size 
    padding = floor(window_size/2); 
    center = padding+1;
    
    right_image = double(padarray(right_image,[padding padding]));
    left_image = double(padarray(left_image,[padding padding]));
    
    % init the weights
    weights = zeros(window_size);
    weights(:) = 1/(window_size * window_size);

    gam_col = 5;
    gam_pos = window_size/2;
    T = 40;
    for row = 1+padding:left_rows+padding
        for col = 1+padding:left_cols+padding
            % get patch from with current pixel as center
            window_left = double(right_image(row-padding:row+padding,col-padding:col+padding));            
            
            % calculate the indexes for sliding window
            max_col = min(col + max_disparity, left_cols + padding);
            min_col = col+min_disparity;
            % init the cost to the worst for each matching algorithm
            best_cost = get_initial_cost(matching_cost);          
                       
            for kk = min_col:max_col
                window_right = double(left_image(row-padding:row+padding,kk-padding:kk+padding));
                
                % compute cost for every pixel of the window
                switch matching_cost
                    case 'SSD'
                        ssd = sum(sum(weights.*double((window_left-window_right).^2)));
                        if ssd < best_cost
                            best_cost = ssd;
                            best_index = kk;
                        end
                    case 'NCC'
                        sum_left = sum(window_left(:).*weights(:));
                        sum_right = sum(window_right(:).*weights(:));

                        sigma_left = sqrt(sum( weights(:).* (window_left(:) - sum_left).^2 ));
                        sigma_right = sqrt(sum( weights(:).* (window_right(:) - sum_right).^2 ));

                        ncc = sum( weights(:).*(window_left(:)-sum_left).*(window_right(:)-sum_right) )/(sigma_left*sigma_right);
                        if ncc > best_cost
                            best_cost = ncc;
                            best_index = kk;
                        end                                         
                       
                    case 'bilateral'
                        num = 0;
                        den = 0;
                        for ll = 1:window_size
                            for mm = 1:window_size
                                p_q = sqrt((ll-padding)^2 + (mm-padding)^2);
                                weight_left = exp(-((abs(window_left(center,center)-window_left(ll,mm))/gam_col) + (p_q/gam_pos)));
                                weight_right = exp(-((abs(window_right(center,center)-window_right(ll,mm))/gam_col) + (p_q/gam_pos)));
                                cost = min(abs(window_left(ll,mm)-window_right(ll,mm)),T);
                                num = num + weight_left*weight_right*cost;
                                den = den + weight_left*weight_right;
                            end
                        end
                        c = num/den;
                        if c < best_cost
                            best_cost = c;
                            best_index = kk;
                        end
                end
            end
            % assign the disparity to the "true" pixel
            disparity(row-padding, col-padding) = best_index-col;
        end        
    end
end


function [val] = get_initial_cost(matching_cost)
    switch matching_cost
        case 'SSD'
            val = Inf;
        case 'NCC'
            val = -Inf;
        case 'bilateral'
            val = Inf;
    end
end