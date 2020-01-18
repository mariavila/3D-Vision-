function [disparity] = stereo_computation(left_image,right_image, ...
    min_disparity, max_disparity, window_size, matching_cost)

[row, col] = size(left_image);
ws_mid = floor(window_size/2);
left_image = padarray(left_image,[ws_mid, ws_mid],'both');
right_image = padarray(right_image,[ws_mid, ws_mid],'both');
weight = ones(window_size)*1/window_size^2;
disparity = zeros(row, col);

for ii = 1+ws_mid:row+ws_mid
    i_min = ii-ws_mid;
    i_max = ii+ws_mid;
    for jj = 1+ws_mid:col+ws_mid
        j_min = jj-ws_mid;
        j_max = jj+ws_mid;
        window_left = left_image(i_min:i_max,j_min:j_max);
        for kk = 1+ws_mid:col+ws_mid
            k_min = kk-ws_mid;
            k_max = kk+ws_mid;
            window_right = right_image(i_min:i_max,k_min:k_max);
            if matching_cost == 'SSD'
                match_cost(kk-ws_mid) = sum(sum(weight.*double((window_left-window_right).^2)));
            end
        end
        [~, d] = sort(match_cost, 'ascend');
        disparity(i_min,j_min) = (d(1)-(j_min));
    end
end

disparity = (disparity-min(disparity(:)));
disparity = disparity/max(disparity(:));
end



