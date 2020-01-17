function [disparity] = stereo_computation(left_image,right_image, ...
    min_disparity, max_disparity, window_size, matching_cost)

[row, col] = size(left_image);
ws_mid = floor(window_size/2);
left_image = padarray(left_image,[ws_mid ws_mid],'both');
right_image = padarray(right_image,[ws_mid ws_mid],'both');
weight = ones(window_size)*1/window_size^2;
disparity = zeros(row, col);

for ii = 1+ws_mid:row+ws_mid
    for jj = 1+ws_mid:col+ws_mid
        window_left = left_image(ii-ws_mid:ii+ws_mid,jj-ws_mid:jj+ws_mid);
        for kk = 1+ws_mid:col+ws_mid
            window_right = right_image(ii-ws_mid:ii+ws_mid,kk-ws_mid:kk+ws_mid);
            if matching_cost == 'SSD'
                match_cost(kk+ws_mid) = sum(sum(weight.*double((window_left-window_right).^2)));
            end
        end
        [~, d] = sort(match_cost, 'descend');
        disparity(ii,jj) = (jj-d(1));
    end
end
end

