function [ output_image ] = apply_H( image, H )
% Applies an homography H to an image 
    %Homogeneous coordinates of the corner points of the image
    size_image = size(image);
    p1 = [1; 1; 1];
    p2 = [1; size_image(2); 1];
    p3 = [size_image(1); 1; 1];
    p4 = [size_image(1); size_image(2); 1];
    
    %Compute location of the corners in the output image
    p1_output_h = H * p1;
    p2_output_h = H * p2;
    p3_output_h = H * p3;
    p4_output_h = H * p4;
    
    %Transform the output corner points to euclidean coordinates
    p1_output_e = [p1_output_h(1)/p1_output_h(3) p1_output_h(2)/p1_output_h(3)];
    p2_output_e = [p2_output_h(1)/p2_output_h(3) p2_output_h(2)/p2_output_h(3)];
    p3_output_e = [p3_output_h(1)/p3_output_h(3) p3_output_h(2)/p3_output_h(3)];
    p4_output_e = [p4_output_h(1)/p4_output_h(3) p4_output_h(2)/p4_output_h(3)];
    
    %Create the output image
    x_min = min([p1_output_e(1) p2_output_e(1) p3_output_e(1) p4_output_e(1)]);
    x_max = max([p1_output_e(1) p2_output_e(1) p3_output_e(1) p4_output_e(1)]);
    y_min = min([p1_output_e(2) p2_output_e(2) p3_output_e(2) p4_output_e(2)]);
    y_max = max([p1_output_e(2) p2_output_e(2) p3_output_e(2) p4_output_e(2)]);
    output_image = zeros(ceil(x_max-x_min), ceil(y_max-y_min), 3);
    
    %Indirect mapping
    [X,Y] = meshgrid(1:size_image(1), 1:size_image(2));
    for i=1:x_max-x_min
        for j=1:y_max-y_min
            p_out_h = [i + x_min; j + y_min; 1];
            p_h = H \ p_out_h;
            p_e = [p_h(1)/p_h(3) p_h(2)/p_h(3)];
            output_image(i,j,1) = interp2(X, Y, image(:,:,1)', [p_e(1)], [p_e(2)]);
            output_image(i,j,2) = interp2(X, Y, image(:,:,2)', [p_e(1)], [p_e(2)]);
            output_image(i,j,3) = interp2(X, Y, image(:,:,3)', [p_e(1)], [p_e(2)]);
        end
    end
end

