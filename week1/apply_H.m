function [ output_image, x_min, x_max, y_min, y_max  ] = apply_H( image, H )
% Applies an homography H to an image 
    %Homogeneous coordinates of the corner points of the image
    size_image = size(image);
    p1 = [1; 1; 1];
    p2 = [1; size_image(1); 1];
    p3 = [size_image(2); 1; 1];
    p4 = [size_image(2); size_image(1); 1];
    
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
    
    %Indirect mapping
    [X,Y] = meshgrid(1:size_image(2), 1:size_image(1));
    [Xo,Yo] = meshgrid(x_min:x_max, y_min:y_max);
    pout = [Xo(:)';Yo(:)';ones(1,numel(Xo))]; 
    pin = inv(H)*pout;
    pin = pin./repmat(pin(3,:),3,1);
  
    r = double(image(:,:,1));
    g = double(image(:,:,2));
    b = double(image(:,:,3));
    output_image = cat(3,...
        interp2(X, Y, r, pin(1,:), pin(2,:), 'cubic', 0),...
        interp2(X, Y, g, pin(1,:), pin(2,:), 'cubic', 0),...
        interp2(X, Y, b, pin(1,:), pin(2,:), 'cubic', 0));
    output_image = uint8(reshape(output_image, [size(Xo,1),size(Xo, 2), 3]));
end
