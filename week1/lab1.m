%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation
%Rotate 45 degrees and translates (10, 10)
H = [0.70710678, 0.70710678, 10;
     -0.70710678, 0.70710678, 10;
     0, 0, 1]; 
I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
H = [1.5, 0.8, 10;
     2, -0.9, 10;
     0, 0, 1];
I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
A = H(1:2, 1:2);
T = H(1:3, 3);
translation = [1 0 T(1);
               0 1 T(2);
               0 0 1];  
[U,S,V] = svd(A);
rotation1 = [U(1,1) U(1,2) 0;
             U(2,1) U(2,2) 0;
             0 0 1];
scale = [S(1,1) S(1,2) 0;
         S(2,1) S(2,2) 0;
         0 0 1];
VT = V';
rotation2 = [VT(1,1) VT(1,2) 0;
             VT(2,1) VT(2,2) 0;
             0 0 1];

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
H2 =  translation * rotation1 * scale * rotation2;
isequal(H, H2)

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I3 = apply_H(I, translation);
I3 = apply_H(I3, rotation1);
I3 = apply_H(I3, scale);
I3 = apply_H(I3, rotation2);
isequal(I2, I3)

%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation
H = [1.5, 0.8, 10;
     -0.3, 0.9, 10;
     0.001, 0.003, 1]; 
I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
%calculate parallel lines
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
%calculate vanishing points
vp1 = cross(l1,l2);
vp1 = vp1/vp1(3);
vp2 = cross(l3,l4);
vp2 = vp2/vp2(3);
%calculate vanishing line
vl = cross(vp1, vp2);

%Write homography normalizing the line
H = [1,0,0;
     0,1,0;
    (vl/vl(3))'];
%H = inv(H);
I2 = apply_H(I, H);
figure; imshow(uint8(I2));
%To transform lines we need the inverse of the transposed homography matrix
H_T = inv(H)';
% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = H_T*l1;
lr2 = H_T*l2;
lr3 = H_T*l3;
lr4 = H_T*l4;
% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
omega = [1 0 0;
         0 1 0;
         0 0 0];
%Angle between lines, we use the scalar product
%Before
disp('Angles of lines before and after affine homography')
disp('Before: ')                     
ang_l12 = radtodeg(acos((omega*l1)'*l2/...
                         (sqrt((omega*l1)'*l1)*sqrt((omega*l2)'*l2)))); 
ang_l13 = radtodeg(acos((omega*l1)'*l3/...
                         (sqrt((omega*l1)'*l1)*sqrt((omega*l3)'*l3)))); 
ang_l24 = radtodeg(acos((omega*l2)'*l4/...
                         (sqrt((omega*l2)'*l2)*sqrt((omega*l4)'*l4))));  
ang_l34 = radtodeg(acos((omega*l3)'*l4/...
                         (sqrt((omega*l3)'*l3)*sqrt((omega*l4)'*l4))));   
fprintf('Angle l1 and l2: %f\n', ang_l12)
fprintf('Angle l1 and l3: %f\n', ang_l13)
fprintf('Angle l2 and l4: %f\n', ang_l24)
fprintf('Angle l3 and l4: %f\n', ang_l34) 

%After                     
disp('After: ')                     
ang_lr12 = radtodeg(acos((omega*lr1)'*lr2/...
                         (sqrt((omega*lr1)'*lr1)*sqrt((omega*lr2)'*lr2)))); 
ang_lr13 = radtodeg(acos((omega*lr1)'*lr3/...
                         (sqrt((omega*lr1)'*lr1)*sqrt((omega*lr3)'*lr3)))); 
ang_lr24 = radtodeg(acos((omega*lr2)'*lr4/...
                         (sqrt((omega*lr2)'*lr2)*sqrt((omega*lr4)'*lr4))));  
ang_lr34 = radtodeg(acos((omega*lr3)'*lr4/...
                         (sqrt((omega*lr3)'*lr3)*sqrt((omega*lr4)'*lr4))));   
fprintf('Angle lr1 and lr2: %f\n', ang_lr12)
fprintf('Angle lr1 and lr3: %f\n', ang_lr13)
fprintf('Angle lr2 and lr4: %f\n', ang_lr24)
fprintf('Angle lr3 and lr4: %f\n', ang_lr34)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

A = [lr1(1)*lr3(1), lr1(1)*lr3(2)+lr1(2)*lr3(1); 
     lr2(1)*lr4(1), lr2(1)*lr4(2)+lr2(2)*lr4(1)];

b = [-lr1(2)*lr3(2), -lr2(2)*lr4(2)];

S = A\b';

S = [S(1), S(2);
     S(2), 1];
K = chol(S);

H = [K(1,1), K(1,1), 0;
     K(2,1), K(2,2), 0;
     0, 0, 1];
I3 = apply_H(I2, inv(H));
figure; imshow(uint8(I3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that are used in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



