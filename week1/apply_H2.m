function [ output_image ] = apply_H2( image, H )
% Applies an homography H to an image 
    %Homogeneous coordinates of the corner points of the image
    size_image = size(image);
%     p =[1 1 1; 1 size_image(1) 1; size_image(2) 1 1; size_image(2) size_image(1) 1]';
    p1 = [1; 1; 1];
    p2 = [1; size_image(1); 1];
    p3 = [size_image(2); 1; 1];
    p4 = [size_image(2); size_image(1); 1];
    
    %Compute location of the corners in the output image
%     po = H*p;
    p1_output_h = H * p1;
    p2_output_h = H * p2;
    p3_output_h = H * p3;
    p4_output_h = H * p4;
    
    %Transform the output corner points to euclidean coordinates
%     po = po./repmat(po(3,:),3,1);
    p1_output_e = [p1_output_h(1)/p1_output_h(3) p1_output_h(2)/p1_output_h(3)];
    p2_output_e = [p2_output_h(1)/p2_output_h(3) p2_output_h(2)/p2_output_h(3)];
    p3_output_e = [p3_output_h(1)/p3_output_h(3) p3_output_h(2)/p3_output_h(3)];
    p4_output_e = [p4_output_h(1)/p4_output_h(3) p4_output_h(2)/p4_output_h(3)];
    
    %Create the output image
    %minim = min(po,[],2)
    %maxim = max(po,[],2)
    x_min = min([p1_output_e(1) p2_output_e(1) p3_output_e(1) p4_output_e(1)]);
    x_max = max([p1_output_e(1) p2_output_e(1) p3_output_e(1) p4_output_e(1)]);
    y_min = min([p1_output_e(2) p2_output_e(2) p3_output_e(2) p4_output_e(2)]);
    y_max = max([p1_output_e(2) p2_output_e(2) p3_output_e(2) p4_output_e(2)]);
    %output_image = zeros(ceil(y_max-y_min), ceil(x_max-x_min), 3);
    
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
    
    
%     for i=1:x_max-x_min
%         for j=1:y_max-y_min
%             p_out_h = [i + x_min; j + y_min; 1];
%             p_h = H \ p_out_h;
%             p_e = [p_h(1)/p_h(3) p_h(2)/p_h(3)];
%             output_image(i,j,1) = interp2(X, Y, image(:,:,1)', [p_e(1)], [p_e(2)]);
%             output_image(i,j,2) = interp2(X, Y, image(:,:,2)', [p_e(1)], [p_e(2)]);
%             output_image(i,j,3) = interp2(X, Y, image(:,:,3)', [p_e(1)], [p_e(2)]);
%         end
%     end


