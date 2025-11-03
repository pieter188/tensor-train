%% Pieter van Klaveren
% 4562461
% this file is used to make the tensor starter kit 1 exercises.
clear all;
close all;
clc
%% Exercise 1, outer product
a = [1;2;3;4];
b = [1;2];
c = [1;2;3];
d = [1;2];
e = [1;2];

% function for calculating the outer product
A1 = outerproduct(a,b,c);
A2 = outerproduct(b,a,c);

% function for transformation of 3-way tensor
A2_t = transformation(A2,a,b,c);
if A2_t == A1
    disp('transformation works correctly')
end
clear A2_t

%% Exercise 2, inner product of two tensors
A = rand(4,2,3,5,6,7);
B = 3*rand(4,2,3,5,6,7);

G = innerproduct(A,B);

%% Exercise 3, mode-n matricization
clear all
d = randi([3,10],1);
A = randi([2,8],[1,d]);
if length(A) ~= d
    disp('error')
end
A = rand(A);

A_mat = Tensor2Mat_moden(A,3); % unfolds the tensor into a matrix

%% Exercise 4, n-mode product
A = rand(5,4,3,2,2,2);
U = diag([1,1,1]);
x=3;
A_n = nmode_product(A,U,x);

%% Exercise 5, tensor-train approximation
clear all; clc;
I = imread('lena.jpg');
% imshow(I);
I = double(I);
% convert matrix into 9-way tensor
I_vec = reshape(I,[numel(I),1]);
I = reshape(I,[4,4,4,4,4,4,4,4,4]);
eps = 0.1;
[B,R] = Tensor2TT_SVD_eps(I,eps);

%% reconstruction of the image via TT
I = double(I);
% convert matrix into 9-way tensor
B_rec = TT2Tensor(B);
B_rec = reshape(B_rec,[512,512]);
B_rec = uint8(B_rec);
imshow(B_rec)

%% Exercise 6, Tuckers method 1 (HOSVD)
% tucker method one, use matricization to compute A_n, which are the left
% eigenvalues of X_n. Use all those to compute G.
I = imread('lena.jpg');
I = double(I);

% convert matrix into 9-way tensor
I_vec = reshape(I,[numel(I),1]);
I = reshape(I,[4,4,4,4,4,4,4,4,4]);
[G,A] = HOSVD(I,R);

%% reconstruction of the image via HOSVD
d = length(size(G));
I_rec_tucker = G;
for a = 1:d
    I_rec_tucker = nmode_product(I_rec_tucker,A{1,a},a);
end
clear a
I_rec_tucker = reshape(I_rec_tucker,[512,512]);
I_rec_tucker = uint8(I_rec_tucker);
imshow(I_rec_tucker)