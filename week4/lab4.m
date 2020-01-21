%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 4: Reconstruction from two views (knowing internal camera parameters) 


addpath('sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Triangulation

% ToDo: create the function triangulate.m that performs a triangulation
%       with the homogeneous algebraic method (DLT)
%
%       The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size

%% Test the triangulate function
% Use this code to validate that the function triangulate works properly

P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X_test = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = euclid(P1 * X_test);
x2_test = euclid(P2 * X_test);

N_test = size(x1_test,2);
X_trian = zeros(4,N_test);
for i = 1:N_test
    X_trian(:,i) = triangulate(x1_test(:,i), x2_test(:,i), P1, P2, [2 2]);
end

% error
euclid(X_test) - euclid(X_trian)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Reconstruction from two views

%% Read images
Irgb{1} = imread('Data/0001_s.png');
Irgb{2} = imread('Data/0002_s.png');
I{1} = sum(double(Irgb{1}), 3) / 3 / 255;
I{2} = sum(double(Irgb{2}), 3) / 3 / 255;
[h,w] = size(I{1});


%% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');


%% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 0.5);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

%vgg_gui_F(Irgb{1}, Irgb{2}, F');

%% Compute candidate camera matrices.

% Camera calibration matrix
K = [2362.12 0 1520.69; 0 2366.12 1006.81; 0 0 1];
scale = 0.3;
H = [scale 0 0; 0 scale 0; 0 0 1];
K = H * K;


% ToDo: Compute the Essential matrix from the Fundamental matrix
E = K'*F*K;

W= [0 -1 0; 
    1 0 0; 
    0 0 1];

[U, D, V] = svd(E);

R1= U*W*V';
if det(R1) < 0
    R1 = -R1;
end
R2= U*W'*V';
if det(R2) < 0
    R2 = -R2;
end

t=U(:,end);

% ToDo: write the camera projection matrix for the first camera
P1 = K*eye(3,4);

% ToDo: write the four possible matrices for the second camera
Pc2 = {};
Pc2{1} = K*[R1 t];
Pc2{2} =  K*[R1 -t];
Pc2{3} =  K*[R2 t];
Pc2{4} = K*[R2 -t];

% HINT: You may get improper rotations; in that case you need to change
%       their sign.
% Let R be a rotation matrix, you may check:
% if det(R) < 0
%     R = -R;
% end

% plot the first camera and the four possible solutions for the second
figure;
plot_camera(P1,w,h);
plot_camera(Pc2{1},w,h);
plot_camera(Pc2{2},w,h);
plot_camera(Pc2{3},w,h);
plot_camera(Pc2{4},w,h);


%% Reconstruct structure
% ToDo: Choose a second camera candidate by triangulating a match.
for i=1:4
    P2 = Pc2{i};
    front = triangulate(x1(:,1), x2(:,1), P1, P2, [w h]);
    proj1 = P1*front;
    proj2 = P2*front;
    if (proj1(3) >= 0) && (proj2(3) >= 0)
        P2 = Pc2{i};
        break
    end
end


% Triangulate all matches.
N = size(x1,2);
X = zeros(4,N);
for i = 1:N
    X(:,i) = triangulate(x1(:,i), x2(:,i), P1, P2, [w h]);
end


%% Plot with colors
r = interp2(double(Irgb{1}(:,:,1)), x1(1,:), x1(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1(1,:), x1(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1(1,:), x1(2,:));
Xe = euclid(X);
figure; hold on;
plot_camera(P1,w,h);
plot_camera(P2,w,h);
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(3,i), -Xe(2,i), 5^2, [r(i) g(i) b(i)]/255, 'filled');
end
axis equal;


%% Compute reprojection error.

% ToDo: compute the reprojection errors
%       plot the histogram of reprojection errors, and
%       plot the mean reprojection error

error1 = gs_errfunction(P1, x1, X);
error2 = gs_errfunction(P2, x2, X);

error1 = error1.^2;
error2 = error2.^2;
histogram([error1 error2])
hold on
mean_err = mean(error1 + error2);
plot([mean_err mean_err], ylim, 'Color','r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Depth map computation with local methods (SSD)

% Data images: 'scene1.row3.col3.ppm','scene1.row3.col4.ppm'
% Disparity ground truth: 'truedisp.row3.col3.pgm'

% Note 1: Use grayscale images
% Note 2: For this first set of images use 0 as minimum disparity 
% and 16 as the the maximum one.

% In this part we ask to implement only the SSD cost
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 21x21,
% 31x31) and the Mean Square Error (MSE). Comment the results.

im_left = imread('Data/scene1.row3.col3.ppm');
im_right = imread('Data/scene1.row3.col4.ppm');
im_left = double(rgb2gray(im_left));
im_right = double(rgb2gray(im_right));


dispSSD_3 = stereo_computation(im_left,im_right, 0, 16, 3, 'SSD');
figure()
imshow(dispSSD_3, [])
%imwrite(uint8(dispSSD_3/max(dispSSD_3(:))*255), 'dispSSD_3.png')

dispSSD_9 = stereo_computation(im_left,im_right, 0, 16, 9, 'SSD');
figure()
imshow(dispSSD_9, [])
%imwrite(uint8(dispSSD_9/max(dispSSD_9(:))*255), 'dispSSD_9.png')

dispSSD_21 = stereo_computation(im_left,im_right, 0, 16, 21, 'SSD');
figure()
imshow(dispSSD_21, [])
%imwrite(uint8(dispSSD_21/max(dispSSD_21(:))*255), 'dispSSD_21.png')

dispSSD_31 = stereo_computation(im_left,im_right, 0, 16, 31, 'SSD');
figure()
imshow(dispSSD_31, [])
%imwrite(uint8(dispSSD_31/max(dispSSD_31(:))*255), 'dispSSD_31.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Depth map computation with local methods (NCC)

% Complete the previous function by adding the implementation of the NCC
% cost.
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 21x21,
% 31x31). Comment the results.

im_left = imread('Data/scene1.row3.col3.ppm');
im_left = double(rgb2gray(im_left));
im_right = imread('Data/scene1.row3.col4.ppm');
im_right = double(rgb2gray(im_right));

dispNCC_3 = stereo_computation(im_left,im_right, 0, 16, 3, 'NCC');
figure()
imshow(dispNCC_3, []);
%imwrite(uint8(dispNCC_3/max(dispNCC_3(:))*255), 'dispNCC_3.png')

dispNCC_9 = stereo_computation(im_left,im_right, 0, 16, 9, 'NCC');
figure()
imshow(dispNCC_9, []);
%imwrite(uint8(dispNCC_9/max(dispNCC_9(:))*255), 'dispNCC_9.png')

dispNCC_21 = stereo_computation(im_left,im_right, 0, 16, 21, 'NCC');
figure()
imshow(dispNCC_21, []);
%imwrite(uint8(dispNCC_21/max(dispNCC_21(:))*255), 'dispNCC_21.png')

dispNCC_31 = stereo_computation(im_left,im_right, 0, 16, 31, 'NCC');
figure()
imshow(dispNCC_31, []);
%imwrite(uint8(dispNCC_31/max(dispNCC_31(:))*255), 'dispNCC_31.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Depth map computation with local methods

% Data images: '0001_rectified_s.png','0002_rectified_s.png'

% Test the functions implemented in the previous section with the facade
% images. Try different matching costs and window sizes and comment the
% results.
% Notice that in this new data the minimum and maximum disparities may
% change.
im_left = imread('Data/0001_rectified_s.png');
im_left = double(rgb2gray(im_left));
im_right = imread('Data/0002_rectified_s.png');
im_right = double(rgb2gray(im_right));

disp_SSD_16 = stereo_computation(im_left,im_right, -50, 20, 3, 'SSD');
figure()
imshow(disp_SSD_16, [])
%imwrite(uint8(disp_SSD_16/max(disp_SSD_16(:))*255), 'dispSSD_16.png')

disp_SSD_50 = stereo_computation(im_left,im_right, 0, 50, 21, 'SSD');
figure()
imshow(disp_SSD_50, [])
%imwrite(uint8(disp_SSD_50/max(disp_SSD_50(:))*255), 'dispSSD_50.png')


disp_SSD_100 = stereo_computation(im_left,im_right, 0, 100, 21, 'SSD');
figure()
imshow(disp_SSD_100, [])
%imwrite(uint8(disp_SSD_100/max(disp_SSD_100(:))*255), 'dispSSD_100.png')

disp_SSD_150 = stereo_computation(im_left,im_right, 0, 150, 21, 'SSD', 'ones');
figure()
imshow(disp_SSD_150, [])
%imwrite(uint8(disp_SSD_150/max(disp_SSD_150(:))*255), 'dispSSD_150.png')


disp_SSD_200 = stereo_computation(im_left,im_right, 0, 200, 21, 'SSD', 'ones');
figure()
imshow(disp_SSD_200, [])
%imwrite(uint8(disp_SSD_200/max(disp_SSD_200(:))*255), 'dispSSD_200.png')


disp_NCC = stereo_computationv2(im_left,im_right, 0, 100, 9, 'NCC', 'ones');
figure()
imshow(disp_NCC, []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Bilateral weights

% Modify the 'stereo_computation' so that you can use bilateral weights (or
% adaptive support weights) in the matching cost of two windows.
% Reference paper: Yoon and Kweon, "Adaptive Support-Weight Approach for Correspondence Search", IEEE PAMI 2006
%
% Comment the results and compare them to the previous results (no weights).
%
% Note: Use grayscale images (the paper uses color images)

im_left = imread('Data/scene1.row3.col3.ppm');
im_left = rgb2gray(im_left);
im_right = imread('Data/scene1.row3.col4.ppm');
im_right = rgb2gray(im_right);

% dispBI_21_30 = stereo_computation(im_left,im_right, 0, 16, 21, 'bilateral');
% figure()
% imshow(dispBI_21_30,[])
% %imwrite(uint8(dispBI_21_30/max(dispBI_21_30(:))*255), 'dispBI_21_30.png')

dispBI_3 = stereo_computation(im_left,im_right, 0, 16, 3, 'bilateral');
figure()
imshow(dispBI_3,[])
%imwrite(uint8(dispBI_3/max(dispBI_3(:))*255), 'dispBI_3.png')

dispBI_9 = stereo_computation(im_left,im_right, 0, 16, 9, 'bilateral');
figure()
imshow(dispBI_9,[])
%imwrite(uint8(dispBI_9/max(dispBI_9(:))*255), 'dispBI_9.png')

dispBI_21 = stereo_computation(im_left,im_right, 0, 16, 21, 'bilateral');
figure()
imshow(dispBI_21,[])
%imwrite(uint8(dispBI_21/max(dispBI_21(:))*255), 'dispBI_21.png')

dispBI_31 = stereo_computation(im_left,im_right, 0, 16, 31, 'bilateral');
figure()
imshow(dispBI_31,[])
%imwrite(uint8(dispBI_31/max(dispBI_31(:))*255), 'dispBI_31.png')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Stereo computation with Belief Propagation

% Use the UGM library used in module 2 and implement a  
% stereo computation method that minimizes a simple stereo energy with 
% belief propagation. 
% For example, use an L2 or L1 pixel-based data term (SSD or SAD) and 
% the same regularization term you used in module 2. 
% Or pick a stereo paper (based on belief propagation) from the literature 
% and implement it. Pick a simple method or just simplify the method they propose.

im_left = imread('Data/scene1.row3.col3.ppm');
im_left = rgb2gray(im_left);
im_right = imread('Data/scene1.row3.col4.ppm');
im_right = rgb2gray(im_right);

% Add  library paths
basedir='UGM/';
addpath(basedir);

%Set model parameters
NumFils = size(im_left,1);
NumCols = size(im_left,2);
K=17; % =number of states of hidden variables
%Pair-wise parameters
smooth_term=[0.0 1.67]; % Potts Model

% Define the unary energy term: data_term
nodePot = unary_potentials(im_left,im_right, 0, 16, 9);

% Create UGM data
[edgePot,edgeStruct] = CreateGridUGMModel(NumFils, NumCols, K ,smooth_term);

% Call UGM inference algorithms
display('Loopy Belief Propagation'); tic;
[nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
[~,c_loopy] = max(nodeBelLBP,[],2);
im_lbp = reshape(c_loopy,size(im_left));toc;

figure
imshow(im_lbp/255, []);xlabel('Loopy Belief Propagation'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth computation with Plane Sweeping
scale_factor = 0.4;

I_scaled{1} = imresize(I{1}, scale_factor);
I_scaled{2} = imresize(I{2}, scale_factor);

P1_scaled = P1;
P2_scaled = P2;
P1_scaled(1:2, :) = P1_scaled(1:2, :)*scale_factor;
P2_scaled(1:2, :) = P2_scaled(1:2, :)*scale_factor;

range_depth = [1 16];
size_window = 21;
cost_function = 'NCC';
threshold = 0.1;
disparity = plane_sweeping(I_scaled, P1, P2, size_window, threshold,cost_function,false);
disparity = imcomplement(disparity);
figure,
imshow(disparity,[])



% Implement the plane sweeping method explained in class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth map fusion 

% In this task you are asked to implement the depth map fusion method
% presented in the following paper:
% B. Curless and M. Levoy. A Volumetric Method for Building Complex
% Models from Range Images. In Proc. SIGGRAPH, 1996.
%
% 1. Use the set of facade images 00xx_s.png to compute depth maps 
% corresponding to different views (and optionally from different pairs of 
% images for the same view).

% 2. Then convert each depth map to a signed distance function defined in 
% a disretized volume (using voxels).
% 3. Average the different signed distance functions, the resulting 
% signed distance is called D.
% 4. Set as occupied voxels (those representing the surface) those 
% where D is very close to zero. The rest of voxels will be considered as 
% empty.
%
% For that you need to compute a depth map from a pair of views in general
% position (non rectified). Thus, you may either use the plane sweep
% algorithm (if you did it) or the local method for estimating depth
% (mandatory task) together with the following rectification method which 
% has an online demo available: 
% http://demo.ipol.im/demo/m_quasi_euclidean_epipolar_rectification/

im_1 = rgb2gray(imread('Data/0001_s.png'));
im_2 = rgb2gray(imread('Data/0002_s.png'));
im_3 = rgb2gray(imread('Data/0003_s.png'));

disp_1_2 = stereo_computationv2(im_1, im_2, 0, 16, 9, 'NCC', 'ones');
disp_2_3 = stereo_computationv2(im_2, im_3, 0, 16, 9, 'NCC', 'ones');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  New view synthesis

% In this task you are asked to implement part of the new view synthesis method
% presented in the following paper:
% S. Seitz, and C. Dyer, View morphing, Proc. ACM SIGGRAPH 1996.

% You will use a pair of rectified stereo images (no need for prewarping
% and postwarping stages) and their corresponding ground truth disparities
% (folder "new_view").
% Remember to take into account occlusions as explained in the lab session.
% Once done you can apply the code to the another pair of rectified images 
% provided in the material and use the estimated disparities with previous methods.

im0 = imread('Data/new_view/im0.png');
disp0 = hdrimread('Data/new_view/disp0.pfm');
im1 = imread('Data/new_view/im1.png');
disp1 = hdrimread('Data/new_view/disp1.pfm');
im0 = rgb2gray(im0);
im1 = rgb2gray(im1);
[rows, cols] = size(disp0);

im0_new = 0*im0;
d0 = 0*disp0;
im1_new = 0*im1;
d1 = 0*disp1;
s = 0.5;

% No disparity
for y=1:rows
    for x=1:cols
        x0 = x-disp0(y,x);
        if (x0>0.5) 
            p0 = round((1-s)*x + s*x0);
            im0_new(y,p0) = (1-s)*im0(y,x) + s*im1(y,round(x0));
            d0(y,p0) = disp0(y,x);
        end      
        x1 = x+disp1(y,x);
        if (x1<cols-0.5)
            p1 = round((1-s)*x1 + s*x);
            im1_new(y,p1) = (1-s)*im1(y,x) + s*im0(y,round(x1));
            d1(y,p1) = disp1(y,x);
        end   
    end
end
figure(3),imshow(im0_new,[])
figure(4),imshow(im1_new,[])
front = (d0-d1)>0;
im_new(front)  = im0_new(front);
im_new(~front) = im1_new(~front);

figure(5),imshow(im_new,[])



% With disparity
im0_new = 0*im0;
im_new = 0*im0;
d0 = 0*disp0;
im1f = fliplr(im1);
disp1f =  fliplr(disp1);
im1_newf = 0*im1;
d1f = 0*disp1;
for y=1:rows
    for x=1:cols
        x0 = x-disp0(y,x);
        if (x0>0.5) 
            p0 = round((1-s)*x + s*x0);
            im0_new(y,p0) = im0(y,x);
            d0(y,p0) = disp0(y,x);
        end      
        
        x1 = x-disp1f(y,x);
        if (x1>0.5)
            p1 = round((1-s)*x + s*x1);
            im1_newf(y,p1) = im1f(y,x);
            d1f(y,p1) = disp1f(y,x);
        end   
    end
end

im1_new =  fliplr(im1_newf);
d1 =  fliplr(d1f);

figure(3),imshow(im0_new,[])
figure(4),imshow(im1_new,[])
front = (d0-d1)>0;
im_new(front)  = im0_new(front);
im_new(~front) = im1_new(~front);

figure(5),imshow(im_new,[])

