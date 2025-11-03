%% Pieter van Klaveren
% 4562461
% this file is used to make the tensor starter kit 2 exercises.
clear all;
close all;
clc;

%% Exercise 1
clear all; clc;
A = rand(4,4,4);
B = rand(4,4,4);
eps = 0.1;
[ATT,Ra] = Tensor2TT_SVD_eps(A,eps);
[BTT,Rb] = Tensor2TT_SVD_eps(A,eps); 

[c,CTT] = innerprodTT(ATT,BTT);

% verify result of Ex1
c_check = sumsqr(A);
if c_check - c < 10^(-10)
    disp('Ex1 done')
end

%% Exercise 2
clear all; clc;
A = rand(3,4,5,2,6,6,3,9);
eps = 0.1;
[ATT,Ra] = Tensor2TT_SVD_eps(A,eps);

k = 5;
ATT_new = sitek(ATT,k);

% verify result of Ex2
normAnew = sumsqr(ATT_new{k,1});
normA = sumsqr(A);

if normAnew - normA < 10^(-10)
    disp('Ex2 done')
end
%% Exercise 3
clear all; clc;
A = rand(3,4,5,6,8,6);
B = rand(3,4,5,6,8,6);
eps = 0.2;
[ATT,Ra] = Tensor2TT_SVD_eps(A,eps);
[BTT,Rb] = Tensor2TT_SVD_eps(A,eps);
clear A B

CTT = addTT(ATT,BTT);

% verify result of Ex3
A_rec = TT2Tensor(ATT);
C_rec = TT2Tensor(CTT);
A_rec = reshape(A_rec,[numel(A_rec),1]);
C_rec = reshape(C_rec,[numel(C_rec),1]);
if norm(C_rec/2 - A_rec) <10^(-10)
    disp('Ex3 done')
end

%% Exercise 4
clear all; clc;
A = rand(3,4,5,6,3,5,6)*10;
eps = 1;
[ATT,Ra] = Tensor2TT_SVD_notrunc(A);

[BTT,Rtt] = roundTT(ATT,eps); 
% denk ook aan computational speed van rounding function.

% verify result of Ex4
A_rec = TT2Tensor(ATT);
B_rec = TT2Tensor(BTT);

if sumsqr(B_rec-A_rec)/sumsqr(A) < eps
    disp('truncation completed succesfully')
end

%% Exercise 5
clear all; clc;
A = rand(72,192);
Rdim = [6,2,2,3];
Cdim = [4,3,8,2];
eps = 10^(-8);

ATT = Mat2TTm_SVD_eps(A,Rdim,Cdim,eps);

A_rec = TTm2Tensor(ATT);
A_rec2 = Tensor2Mat(A_rec,4,4);

% % compute the transpose of A via TTm
% for i = 1:length(dim)/2
%     ATT{i,1} = permute(ATT{i,1},[1,3,2,4]);
% end



%% Exercise 6, matrix-vector product
clear all; clc;
A = rand(1920,500);         % create A matrix
x = rand(500,1);            % create x column vector
B = A*x;                    % compute b (column vector)
RdimA = [4,8,10,6];         % row dimension of A
CdimA = [10,5,5,2];         % column dimension of A
dimA = [RdimA,CdimA];       % TTm: [up, under], the desired 
                            %  dimensions, n(k), of the TTm-cores of A.
dimx = [10,5,5,2];          % TT: [up], the desired dimensions, n(k), of
                            %  the TT-cores of x.                          
eps = 10^(-7);
x_tensor = reshape(x,dimx); % reshape the vector x into a tensor, according
                            %  to the dimensions: dimx, so that it can be
                            %  used in the TTSVD algorithm.
ATT = Mat2TTm_SVD_eps(A,RdimA,CdimA,eps);       % create the TTm of A.
xTT = Tensor2TT_SVD_eps(x_tensor,eps);      % create the TT of x. Checked.

BTT = matrixvec(ATT,xTT);

% Check if the reconstructed vector B_rec is equal to the original vector B
B_rec = TT2Tensor(BTT);
B_rec = reshape(B_rec,[numel(B_rec),1]);

