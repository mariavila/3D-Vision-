function [disparity] = stereo_computation(left_image,right_image, min_disparity, max_disparity, window_size, matching_cost, weight_f)

    disparity = zeros(size(left_image));
    [left_rows, left_cols] = size(left_image);
    
    % we need to pad half the window size 
    padding = floor(window_size/2); 
    
    left_image = padarray(left_image,[padding padding]);
    right_image = padarray(right_image,[padding padding]);
    
    % init the weights
    switch weight_f
    case 'ones'
        weights = zeros(window_size);
        weights(:) = 1/(window_size * window_size);
    case 'gaussian'
        x = -padding:padding;
        [x,y] = meshgrid(x,x);
        sigma = 3;
        exponent = ((x).^2 + (y).^2)./(2*sigma^2);
        amplitude = 1 / (sigma * sqrt(2*pi));  
        weights = amplitude  * exp(-exponent);
    end


    gam_col = 7;
    gam_pos = 7;
    for row = 1+padding:left_rows+padding
        for col = 1+padding:left_cols+padding
            % get patch from with current pixel as center
            window_left = double(left_image(row-padding:row+padding,col-padding:col+padding));            
            
            % calculate the indexes for sliding window
            min_col = max(max(1+padding, min_disparity), col - max_disparity);
            max_col = min(left_cols + padding, col + max_disparity);
            
            % init the cost to the worst for each matching algorithm
            best_cost = get_initial_cost(matching_cost);
            
            if strcmp(weight_f, 'bilateral')
                d_c = abs(ones(window_size)*window_left(padding,padding) - window_left);
                qq = -padding:padding;
                [qq1,qq2] = meshgrid(qq,qq);                
                d_g = sqrt((0-qq1)^2+(0-qq2^2));
                weights = exp(-d_c/gam_col).*exp(-d_g/gam_pos);
            end
            
            
            for kk = min_col:max_col
                window_right = double(right_image(row-padding:row+padding,kk-padding:kk+padding));
                
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
                end
            end
            % assign the disparity to the "true" pixel
            disparity(row-padding, col-padding) = abs(col-best_index);
        end        
    end
end


function [val] = get_initial_cost(matching_cost)
    switch matching_cost
        case 'SSD'
            val = Inf;
        case 'NCC'
            val = -Inf;     
    end
end