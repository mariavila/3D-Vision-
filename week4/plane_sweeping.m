function [ I_depth ] = plane_sweeping( I,P1,P2,window_size, matching_function )

weights = ones(window_size)/window_size.^2;
[rows,cols] = size(I{1});
corners = [1, cols, 1, rows];
I_depth = inf(rows,cols);
best_matching = inf(rows,cols);
if strcmpi(matching_function,'NCC')
    best_matching=-best_matching;
end
pad = floor(window_size/2);
I1 = padarray(I{1},[pad,pad],'replicate');
for depth=1:20
    fronto_parallel_plane = (P2(3,:)-[0,0,0,depth])';
    
    A = inv([P2; fronto_parallel_plane']);
    H = P1*A(:,1:3);
    
    I2warped = apply_H_v2(I{2}, H, corners);
    
    I2warped = padarray(I2warped,[pad,pad],'replicate');
    I2warped(I2warped==0)=Inf;
    matching = zeros(rows,cols);
    I1 = permute(I1,[1,3,2]);
    I2warped = permute(I2warped,[1,3,2]);
    for i=1+pad:rows+pad
        I1_i = zeros(window_size,window_size,cols);
        I2_i = zeros(window_size,window_size,cols);
        for wj=1:window_size
            I1_i(:,wj,:) = I1(i-pad:i+pad,wj:cols+(2*pad)-window_size+wj);
            I2_i(:,wj,:) = I2warped(i-pad:i+pad,wj:cols+(2*pad)-window_size+wj);
        end
        if strcmpi(matching_function,'SSD')
            matching(i-pad,:) = permute(sum(sum(weights.*(I1_i-I2_i).^2,1),2),[1,3,2]);
        elseif strcmpi(matching_function,'NCC')
            sum_I1_row = sum(sum(weights .* I1_i,1),2);
            sum_I2_row = sum(sum(weights .* I2_i,1),2);
            
            sigma_1 = sqrt(sum(sum(weights .* (I1_i - sum_I1_row).^2,1),2));
            sigma_2 = sqrt(sum(sum(weights .* (I2_i - sum_I2_row).^2,1),2));
            
           
            aux = sum(sum(weights .* (I1_i - sum_I1_row) .* (I2_i - sum_I2_row), 1), 2) ./ (sigma_1.*sigma_2);
            matching(i-pad,:)=aux;
        end
    end
    
    if strcmpi(matching_function,'SSD')
        I_depth(best_matching > matching) = depth;
        best_matching(best_matching > matching) = matching(best_matching > matching);
    elseif strcmpi(matching_function,'NCC')
        I_depth(best_matching < matching) = depth;
        best_matching(best_matching < matching) = matching(best_matching < matching);
    end
end
end