%% Concatenating matrices in TT-format
% C = [A 0] + [0 B]
clear all; clc

%% [A 0]
A = rand(20);
V = [1 0];
C = kron(V,A);

% partitioning of row and column of A
Row = [5 2 2];
Col = [2 5 2];
% tensor A according to row.*col
Aten = Mat2Tensor(A,Row,Col);
eps = 0.01;
ATT = Tensor2TT_SVD_eps(Aten,eps);
ATTm = TT2TTm(ATT,Row,Col);
clear A Aten 

%% connect V to the last core of A
a = [size(ATTm{3,1}),1];
ATTm{3,1} = reshape(ATTm{3,1},[a(1)*a(2)*a(3),a(4)]);
ATTm{3,1} = ATTm{3,1}*V;
ATTm{3,1} = reshape(ATTm{3,1},[a(1),a(2),a(3),2]);
ATTm{3,1} = reshape(ATTm{3,1},[a(1),a(2),a(3)*2,a(4)]);

%% check if Cmat equals C
[CTT,Rdim,Cdim] = TTm2TT(ATTm);
Cten = TT2Tensor(CTT);
Cmat = Tensor2Mat(Cten,Rdim,Cdim);
if norm(C-Cmat)<10^(-10)
   display('correct') 
end

%% [0 B]

B = rand(20);
W = [0 1];
D = kron(W,B);

% partitioning of row and column of A
Row = [5 2 2];
Col = [2 5 2];
% tensor A according to row.*col
Bten = Mat2Tensor(B,Row,Col);
eps = 0.01;
BTT = Tensor2TT_SVD_eps(Bten,eps);
BTTm = TT2TTm(BTT,Row,Col);
clear B Bten 

%% connect V to the last core of A
a = [size(BTTm{3,1}),1];
BTTm{3,1} = reshape(BTTm{3,1},[a(1)*a(2)*a(3),a(4)]);
BTTm{3,1} = BTTm{3,1}*W;
BTTm{3,1} = reshape(BTTm{3,1},[a(1),a(2),a(3),2]);
BTTm{3,1} = reshape(BTTm{3,1},[a(1),a(2),a(3)*2,a(4)]);

%% check if Cmat equals C
[DTT,Rdim,Cdim] = TTm2TT(BTTm);
Dten = TT2Tensor(DTT);
Dmat = Tensor2Mat(Dten,Rdim,Cdim);
if norm(D-Dmat)<10^(-10)
   display('correct') 
end

%% addition of ATT and BTT
[ATT,rA,cA] = TTm2TT(ATTm);
[BTT,rB,cB] = TTm2TT(BTTm);

CTT = addTT(ATT,BTT);

%% check whether CTTm equals C+D
C_check = Tensor2Mat(TT2Tensor(TTm2TT(CTTm)),rA,cA);
if norm(C_check-(C+D))<10^(-10)
   display('correct') 
end

%% check for [L Q]
clear all
clc

n = 9216;
% partitioning of row and column of A
Row = [4 2 4 2 4 3 4 3];
Col = [4 2 4 2 4 3 4 3];

% initial square-root covariance matrix
P = eye(n)*0.1;
L = chol(P);
clear P
% from matrix L to TTm
Lten = Mat2Tensor(L,Row,Col);
eps = 0.01;
LTT = Tensor2TT_SVD_eps(Lten,eps);
LTTm = TT2TTm(LTT,Row,Col);
clear L Lten LTT

% process noise covariance matrix
Q = chol(eye(n)*5);
Qten = Mat2Tensor(Q,Row,Col);
eps = 0.01;
QTT = Tensor2TT_SVD_eps(Qten,eps);
QTTm = TT2TTm(QTT,Row,Col);
clear Q Qten QTT

%% concatTT
if size(LTTm) ~= size(QTTm)
    error('The matrix dimensions of A and B are not the same')
end 

%% [LTTm 0]
[d,~] = size(LTTm);
V = [1 0];
% connect V to the last core of A
a = [size(LTTm{d,1}),1];
LTTm{d,1} = reshape(LTTm{d,1},[a(1)*a(2)*a(3),a(4)]);
LTTm{d,1} = LTTm{d,1}*V;
LTTm{d,1} = reshape(LTTm{d,1},[a(1),a(2),a(3),2]);
LTTm{d,1} = reshape(LTTm{d,1},[a(1),a(2),a(3)*2,a(4)]);

% TTm -> TT
[LTT,nrL,ncL] = TTm2TT(LTTm);

%% [0 QTTm]
W = [0 1];
% connect V to the last core of A
a = [size(QTTm{d,1}),1];
QTTm{d,1} = reshape(QTTm{d,1},[a(1)*a(2)*a(3),a(4)]);
QTTm{d,1} = QTTm{d,1}*W;
QTTm{d,1} = reshape(QTTm{d,1},[a(1),a(2),a(3),2]);
QTTm{d,1} = reshape(QTTm{d,1},[a(1),a(2),a(3)*2,a(4)]);

% TTm -> TT
[QTT,nrQ,ncQ] = TTm2TT(QTTm);

%% [A 0] + [0 B]
% check whether the dimensions are the same
if ~isequal(nrL,nrQ) || ~isequal(ncL,ncQ)
   error('The dimensions the two TT are not identical') 
end
% addTT is done with TT and not with TTm
CTT = addTT(LTT,QTT);
CTTm = TT2TTm(CTT,nrL,ncL);
