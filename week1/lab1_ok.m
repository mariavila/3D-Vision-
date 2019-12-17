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
% a = [1 0 0;
%      0 -1 0;
%      0 0 0];
% I2 = apply_H(I, a*H*a);
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

i = 227;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 534;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
%calculate parallel lines
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);


c1 = cross(l1, l3);
c1 = c1/c1(3);
c2 = cross(l1, l4);
c2 = c2/c2(3);
c3 = cross(l2, l3);
c3 = c3/c3(3);
c4 = cross(l2, l4);
c4 = c4/c4(3);

l5 = cross(c1,c4);
l6 = cross(c2,c3);


% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(c1(1), c1(2)); 
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'y');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
% First we calculate two vanishing points from 2 pairs of parallel lines
% and normalize them, make the third coordinate = 0.
vp1 = cross(l1,l2);
vp1 = vp1/vp1(3);
vp2 = cross(l3,l4);
vp2 = vp2/vp2(3);

% From the vanishing points calculate vanishing line vl
vl = cross(vp1, vp2);
% This line is invariant under affine transformations.
H = [1,0,0;
     0,1,0;
    (vl/vl(3))'];
I2 = apply_H(I, H);
%temp = maketform('projective',H');
%I2 = imtransform(I,temp);
%figure; imshow(uint8(I2));
figure; imshow(uint8(I2));
H_T = inv(H)';


lr1 = H_T*l1;
lr2 = H_T*l2;
lr3 = H_T*l3;
lr4 = H_T*l4;
lr5 = H_T*l5;
lr6 = H_T*l6;
omega = [1 0 0;
         0 1 0;
         0 0 0];
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


figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'y');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'y');


A = [lr1(1)*lr3(1), (lr1(1)*lr3(2)+lr1(2)*lr3(1));
     lr5(1)*lr6(1), (lr5(1)*lr6(2)+lr5(2)*lr6(1))];
b = [-lr1(2)*lr3(2); -lr5(2)*lr6(2)];

S = A\b;

S = [S(1), S(2);
    S(2), 1];

K = chol(S);
K = inv(K');
H = [K(1,1), K(1,2), 0;
     K(2,1), K(2,2), 0;
     0, 0, 1];
I3 = apply_H(I2, H);
H_T = inv(H)';

lr1 = H_T*lr1;
lr2 = H_T*lr2;
lr3 = H_T*lr3;
lr4 = H_T*lr4;
lr5 = H_T*lr5;
lr6 = H_T*lr6;
omega = [1 0 0;
         0 1 0;
         0 0 0];
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

figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'y');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'y');


I=imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');
I = I(:,1:size(I,2)/2,:);
% indices of lines
i = 159;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 614;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 541;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 645;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';


% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
%calculate parallel lines
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

c1 = cross(l1, l3);
c1 = c1/c1(3);
c2 = cross(l1, l4);
c2 = c2/c2(3);
c3 = cross(l2, l3);
c3 = c3/c3(3);
c4 = cross(l2, l4);
c4 = c4/c4(3);

l5 = cross(c1,c4);
l6 = cross(c2,c3);


% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(c1(1), c1(2)); 
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'y');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
% First we calculate two vanishing points from 2 pairs of parallel lines
% and normalize them, make the third coordinate = 0.
vp1 = cross(l1,l2);
vp1 = vp1/vp1(3);
vp2 = cross(l3,l4);
vp2 = vp2/vp2(3);

% From the vanishing points calculate vanishing line vl
vl = cross(vp1, vp2);
% This line is invariant under affine transformations.
H = [1,0,0;
     0,1,0;
    (vl/vl(3))'];
[I2,x_min0, x_max0, y_min0, y_max0 ] = apply_H(I, H);
%temp = maketform('projective',H');
%I2 = imtransform(I,temp);
%figure; imshow(uint8(I2));
figure; imshow(uint8(I2));
H_T = inv(H)';


lr1 = H_T*l1;
lr2 = H_T*l2;
lr3 = H_T*l3;
lr4 = H_T*l4;
lr5 = H_T*l5;
lr6 = H_T*l6;

figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'y');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'y');

A = [lr1(1)*lr3(1), (lr1(1)*lr3(2)+lr1(2)*lr3(1));
     lr5(1)*lr6(1), (lr5(1)*lr6(2)+lr5(2)*lr6(1))];
b = [-lr1(2)*lr3(2); -lr5(2)*lr6(2)];

S = A\b;

S = [S(1), S(2);
    S(2), 1];

K = chol(S);
K = inv(K');
%K = K';
H = [K(1,1), K(1,2), 0;
     K(2,1), K(2,2), 0;
     0, 0, 1];
[I3,x_min, x_max, y_min, y_max ] = apply_H(I2, H);
H(1,3) = -x_min;
H(2,3) = -y_min;
H_T = inv(H)';
lrr1 = H_T*lr1;
lrr2 = H_T*lr2;
lrr3 = H_T*lr3;
lrr4 = H_T*lr4;
lrr5 = H_T*lr5;
lrr6 = H_T*lr6;

omega = [1 0 0;
         0 1 0;
         0 0 0];
%After                     
disp('After: ')                     
ang_lr12 = radtodeg(acos((omega*lrr1)'*lrr2/...
                         (sqrt((omega*lrr1)'*lrr1)*sqrt((omega*lrr2)'*lrr2)))); 
ang_lr13 = radtodeg(acos((omega*lrr1)'*lrr3/...
                         (sqrt((omega*lrr1)'*lrr1)*sqrt((omega*lrr3)'*lrr3)))); 
ang_lr24 = radtodeg(acos((omega*lrr2)'*lrr4/...
                         (sqrt((omega*lrr2)'*lrr2)*sqrt((omega*lrr4)'*lrr4))));  
ang_lr34 = radtodeg(acos((omega*lrr3)'*lrr4/...
                         (sqrt((omega*lrr3)'*lrr3)*sqrt((omega*lrr4)'*lr4))));   
fprintf('Angle lr1 and lr2: %f\n', ang_lr12)
fprintf('Angle lr1 and lr3: %f\n', ang_lr13)
fprintf('Angle lr2 and lr4: %f\n', ang_lr24)
fprintf('Angle lr3 and lr4: %f\n', ang_lr34)


figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t, -(lrr2(1)*t + lrr2(3)) / lrr2(2), 'y');
plot(t, -(lrr3(1)*t + lrr3(3)) / lrr3(2), 'y');
plot(t, -(lrr4(1)*t + lrr4(3)) / lrr4(2), 'y');
plot(t, -(lrr5(1)*t + lrr5(3)) / lrr5(2), 'y');
plot(t, -(lrr6(1)*t + lrr6(3)) / lrr6(2), 'y');
