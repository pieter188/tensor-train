%% Concatenating matrices in TT-format and QR decomposition.
% Perform QR decomposition via (Modified) Gram-Schmidt procedure. Find the 
clear all; clc
%% generate matrices
L = tril(rand(10));
row = [5 2];
col = [5 2];
Lten = Mat2Tensor(L,row,col);
eps = 0.01;
LTT = Tensor2TT_SVD_eps(Lten,eps);
LTTm = TT2TTm(LTT,row,col);
clear Lten LTT

C = [0 0 1 0 0 0 0 0 0 0];
rowC = [1 1];
colC = [5 2];
Cten = Mat2Tensor(C,rowC,colC);
eps = 0.01;
CTT = Tensor2TT_SVD_eps(Cten,eps);
CTTm = TT2TTm(CTT,rowC,colC);
clear Cten CTT

W  = [C*L;L];
Wcom = [zeros(9,10);C*L;L];

CL = MatMatTT(CTTm,LTTm);

%% check if CL is correct -> yes it is!
% [CL,rowC,colC] = TTm2TT(CL);
% CL = TT2Tensor(CL);
% CL = Tensor2Mat(CL,rowC,colC);

%% 
