%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

addpath('sift');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images

imargb = imread('Data/llanes/llanes_a.jpg');
imbrgb = imread('Data/llanes/llanes_b.jpg');
imcrgb = imread('Data/llanes/llanes_c.jpg');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

figure;
imshow(imargb);%image(imargb)
hold on;
plot(points_a(1,:), points_a(2,:),'+y');
figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(points_b(1,:), points_b(2,:),'+y');
figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
    matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
    matches_bc(:,inliers_bc), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

corners = [-400 1200 -100 650];
Hbb = eye(3);
iwb = apply_H_v2(imbrgb, Hbb , corners);   % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab, corners);    % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);    % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');
% DISCUSSION: In this case the method to build the panorama works 
% because we can assume the camera that took the images is static 
% as the scene is far from the camera. For each image the camera 
% only changes the rotation.

% ToDo: compute the mosaic with castle_int images
imargb = imread('Data/castle_int/0016_s.png');
imbrgb = imread('Data/castle_int/0015_s.png');
imcrgb = imread('Data/castle_int/0014_s.png');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

matches_ab = siftmatch(desc_a, desc_b);
matches_bc = siftmatch(desc_b, desc_c);

th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

corners = [-570 1550 -250 780];
Hbb = eye(3);
iwb = apply_H_v2(imbrgb, Hbb , corners);   
iwa = apply_H_v2(imargb, Hab, corners);    
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);   

figure;
imshow(max(iwc, max(iwb, iwa)));
title('Mosaic A-B-C');
% DISCUSSION: in this case the tractor that appears on images ABC
% is impossible to reconstruct after the transformations. It is 
% because it is very close to the camera and we have occlusions 
% which makes it impossible to be matched. The facade on the other
% hand can be built into a mosaic as it fulfills the conditions.

% ToDo: compute the mosaic with aerial images set 13
imargb = imread('Data/aerial/site13/frame00000.png');
imbrgb = imread('Data/aerial/site13/frame00002.png');
imcrgb = imread('Data/aerial/site13/frame00003.png');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

matches_ab = siftmatch(desc_a, desc_b);
matches_bc = siftmatch(desc_b, desc_c);

th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

corners = [-300 1300 -100 950];
Hbb = eye(3);
iwb = apply_H_v2(imbrgb, Hbb , corners);   
iwa = apply_H_v2(imargb, Hab, corners);    
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);   

figure;
imshow(max(iwc, max(iwb, iwa)));
title('Mosaic A-B-C');

% DISCUSSION: see that the inferior left corner of the result 
% image is a little sloppy. We have studied why it might be 
% and for some reason the ransac discards the matches in that 
% zone so there the homography is not exact and the result image
% gives this double building effect in that zone. We think it might 
% be something similar to the case of the tractor but less severe. 
% This case might be in the limit of a good performance of the method.
% See the document to see the images that ilustrate the discussion.

% ToDo: compute the mosaic with aerial images set 22
imargb = double(imread('Data/aerial/site22/frame_00001.tif'));
imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'));
imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'));

ima = imargb;
imb = imbrgb;
imc = imcrgb;

[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

matches_ab = siftmatch(desc_a, desc_b);
matches_bc = siftmatch(desc_b, desc_c);

th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

corners = [-400 1450 -100 1050];
Hbb = eye(3);
iwb = apply_H_v2(imbrgb, Hbb , corners);   
iwa = apply_H_v2(imargb, Hab, corners);    
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);   

figure;
imshow(max(iwc, max(iwb, iwa)));
title('Mosaic A-B-C');

% DISCUSSION: between images the camera has clearly changed the position 
% and it is the main hypothesis in order for this method to work. That
% is why the result cannot build the mosaic properly. For the objects
% closer to the ground (streets, roads, the river,...) the camera is very 
% far from them and we can assume we have an static camera and they are 
% correctly transformed so they match. Also mention that it takes quite a 
% long time to calculate this mosaic and we think that is because of the 
% lack of correct matches and it takes longer for the RANSAC to compute.
% See the document to see the images that ilustrate the discussion.

% ToDo: comment the results in every of the four cases: hypothetise why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm

% Homography ab

x = points_a(1:2, matches_ab(1,inliers_ab));  %ToDo: set the non-homogeneous point coordinates of the 
xp = points_b(1:2, matches_ab(2,inliers_ab)); %      point correspondences we will refine with the geometric method
% x = [points_a(1:2, matches_ab(1,inliers_ab)); ones(1, length(matches_ab(1,inliers_ab)))];
% xp = [points_b(1:2, matches_ab(2,inliers_ab)); ones(1, length(matches_ab(1,inliers_ab)))];
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction_v2( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction_v2(t, Xobs), P0, [], [], options);

Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

vgg_gui_H(imargb, imbrgb, Hab_r);

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = P(10:end);
xhat = reshape(xhat, [2,size(xhat,1)/2]); % from nx1 to 2x(n/2)
xhat = [xhat ; ones(1,size(xhat,2))]; % from euclidean to homogeneous
xhatp = Hab_r*xhat;

xhat = euclid(xhat);
xhatp = euclid(xhatp);

figure;
imshow(imargb, []);%image(imargb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure;
imshow(imbrgb, []);%image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');

%%  Homography bc

% ToDo: refine the homography bc with the Gold Standard algorithm

x = points_b(1:2, matches_bc(1,inliers_bc));  %ToDo: set the non-homogeneous point coordinates of the 
xp = points_c(1:2, matches_bc(2,inliers_bc)); %      point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [Hbc(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction_v2( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction_v2(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction_v2( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

vgg_gui_H(imbrgb, imcrgb, Hbc_r);


%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = P(10:end);
xhat = reshape(xhat, [2,size(xhat,1)/2]); % from nx1 to 2x(n/2)
xhat = [xhat ; ones(1,size(xhat,2))]; % from euclidean to homogeneous
xhatp = Hbc_r*xhat;

xhat = euclid(xhat);
xhatp = euclid(xhatp);

figure;
imshow(imbrgb, []);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure;
imshow(imcrgb, []);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');

%% Build mosaic
corners = [-400 1450 -100 1050];
iwb = apply_H_v2(imbrgb, eye(3), corners); % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab_r, corners); % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc_r), corners); % ToDo: complete the call to the function
% iwc = zeros(size(iwc), 'uint8');
figure;
imshow(max(iwc, max(iwb, iwa)), []);%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern

clear all;

%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = pointsT(1:2, matches(1, :));
    x1 = [x1;ones(1,size(x1,2))];
    x2 = points{i}(1:2, matches(2, :));
    x2 = [x2;ones(1,size(x2,2))];
    H{i} = 0;
%     [H{i}, inliers] =  ransac_homography_adaptive_loop(homog(x1), homog(x2), 3, 1000);
    [H{i}, inliers] =  ransac_homography_adaptive_loop(x1, x2, 3, 1000);
    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

%     Play with the homography
%     vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic
n_ims = size(I,2);
vij = zeros(2*n_ims,6);
count = 1;
for i = 1:n_ims
    %v_12
    A(count,:) = [H{i}(1,1)*H{i}(1,2), H{i}(1,1)*H{i}(2,2)+H{i}(2,1)*H{i}(1,2),...
                  H{i}(1,1)*H{i}(3,2)+H{i}(3,1)*H{i}(1,2), H{i}(2,1)*H{i}(2,2),...
                  H{i}(2,1)*H{i}(3,2)+H{i}(3,1)*H{i}(2,2), H{i}(3,1)*H{i}(3,2)];
    %v_11
    v_11 = [H{i}(1,1)*H{i}(1,1), H{i}(1,1)*H{i}(2,1)+H{i}(2,1)*H{i}(1,1),...
                  H{i}(1,1)*H{i}(3,1)+H{i}(3,1)*H{i}(1,1), H{i}(2,1)*H{i}(2,1),...
                  H{i}(2,1)*H{i}(3,1)+H{i}(3,1)*H{i}(2,1), H{i}(3,1)*H{i}(3,1)];
    %v_22
    v_22 = [H{i}(1,2)*H{i}(1,2), H{i}(1,2)*H{i}(2,2)+H{i}(2,2)*H{i}(1,2),...
                  H{i}(1,2)*H{i}(3,2)+H{i}(3,2)*H{i}(1,2), H{i}(2,2)*H{i}(2,2),...
                  H{i}(2,2)*H{i}(3,2)+H{i}(3,2)*H{i}(2,2), H{i}(3,2)*H{i}(3,2)];
              
    %v_11-v_22
    A(count+1,:) = v_11-v_22;
    count = count + 2;
end
[~,~,V] = svd(A);
X = V(:,end);
w = [X(1), X(2), X(3); X(1), X(4), X(5); X(3), X(5), X(6)]; % ToDo
 
%% Recover the camera calibration.

K = inv(chol(w)); % ToDo
% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    % ToDo: compute r1, r2, and t{i}
    %r1 = ...
    %r2 = ...
    %t{i} = ...
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U S V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
% figure; hold;
% plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
% for i = 1:N
%     vgg_scatter_plot( [...   ...   ...   ...   ...], 'r');
% end

%% Augmented reality: Plot some 3D points on every camera.
%[Th, Tw] = size(Tg);
%cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

%X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));

%for i = 1:N
%     figure; colormap(gray);
%     imagesc(Ig{i});
%     hold on;
%     x = euclid(P{i} * homog(X));
%     vgg_scatter_plot(x, 'g');
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
%%              DLT algorithm (folder "logos").
%%              Interpret and comment the results.
imrgb_building = imread('Data/logos/UPFbuilding.jpg');
imrgb_stand = imread('Data/logos/UPFstand.jpg');
imrgb_logo_upf = imread('Data/logos/logoUPF.png');

im_building = sum(double(imrgb_building), 3) / 3 / 255;
im_stand = sum(double(imrgb_stand), 3) / 3 / 255;
im_logo_upf = sum(double(imrgb_logo_upf), 3) / 3 / 255;

%Compute SIFT keypoints
[points_building, desc_building] = sift(im_building, 'Threshold', 0.01);
[points_stand, desc_stand] = sift(im_stand, 'Threshold', 0.01);
[points_logo_upf, desc_logo_upf] = sift(im_logo_upf, 'Threshold', 0.01);

%% Match SIFT keypoints 
% between logo and building
matches_logo_building = siftmatch(desc_logo_upf, desc_building);

% between logo and stand
matches_logo_stand = siftmatch(desc_logo_upf, desc_stand);

%%Compute the homography using the DLT algorithm
%Compute homography (normalized DLT) between logo and building
th = 3;
xlb_logo = [points_logo_upf(1:2, matches_logo_building(1,:)); ones(1, length(matches_logo_building))];
xlb_building = [points_building(1:2, matches_logo_building(2,:)); ones(1, length(matches_logo_building))];
[H_logo_building, inliers_lb] = ransac_homography_adaptive_loop(xlb_logo, xlb_building, th, 1000);

%Compute homography (normalized DLT) between logo and stand
th = 3;
xls_logo = [points_logo_upf(1:2, matches_logo_stand(1,:)); ones(1, length(matches_logo_stand))];
xls_stand = [points_stand(1:2, matches_logo_stand(2,:)); ones(1, length(matches_logo_stand))];
[H_logo_stand, inliers_ls] = ransac_homography_adaptive_loop(xls_logo, xls_stand, th, 1000);


%%Show logo location
[sy, sx, ~] = size(im_logo_upf);
corners = [0 sx sx 0; 0 0 sy sy];
x = [corners ; ones(1,4)];

%In building image 
corners_building = H_logo_building * x;
corners_building = corners_building(1:2, :) ./ repmat(corners_building(3,:),2,1);
figure;
imshow(imrgb_building);
hold on;
line([corners_building(1,:) corners_building(1,1)], [corners_building(2,:) corners_building(2,1)], 'color', 'y', 'LineWidth', 2);

%In stand image 
corners_stand = H_logo_stand * x;
corners_stand = corners_stand(1:2, :) ./ repmat(corners_stand(3,:),2,1);
figure;
imshow(imrgb_stand);
hold on;
line([corners_stand(1,:) corners_stand(1,1)], [corners_stand(2,:) corners_stand(2,1)], 'color', 'y', 'LineWidth', 2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Replace the logo of the UPF by the master logo
%%              in one of the previous images using the DLT algorithm.
imrgb_logo_master = imread('Data/logos/logo_master.png');
im_logo_master = sum(double(imrgb_logo_master), 3) / 3 / 255;

[sy, sx, ~] = size(im_logo_master);
[sy2, sx2, ~] = size(im_logo_upf);

scale = [sx2/sx 0 0; 0 sy2/sy 0; 0 0 1];
H = H_logo_stand * scale;

corners = [0 sx 0 sy];
imt_master = apply_H_v2(imrgb_logo_master, H, corners); 
Hss = eye(3);  
imt_stand = apply_H_v2(imrgb_stand, Hss, corners);
figure;
imshow(max(imt_stand, imt_master));
%imshow(imt_stand);
%hold on;
%imshow(imt_master);




