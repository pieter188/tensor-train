function C = concatTTm(A,B)
% In the first section this function creates two structures [A 0] and [0
% B]. It is required that A and B have the same dimension in order to
% succesfully do the concatenation.
% INPUT: A and B, two TTm's
% OUTPUT: C, TT concatenation of A and B

if size(A.dim,1) ~= size(B.dim,1)
    error('The amount of cores of A and B are not the same')
end 

%% [A 0]
d = size(A.dim,1);
V = [1 0];
% connect V to the last core of A
a = A.dim(d,:);
A.cores{d,1} = reshape(A.cores{d,1},[a(1)*a(2)*a(3),a(4)]);
A.cores{d,1} = A.cores{d,1}*V;
A.cores{d,1} = reshape(A.cores{d,1},[a(1),a(2),a(3),2]);
A.cores{d,1} = reshape(A.cores{d,1},[a(1),a(2),a(3)*2,a(4)]);
A.dim(d,:) = [a(1),a(2),a(3)*2,a(4)];

% TTm -> TT
[ATT,nrA,ncA] = TTm2TT(A);

%% [0 B]
W = [0 1];
% connect V to the last core of A
a = B.dim(d,:);
B.cores{d,1} = reshape(B.cores{d,1},[a(1)*a(2)*a(3),a(4)]);
B.cores{d,1} = B.cores{d,1}*W;
B.cores{d,1} = reshape(B.cores{d,1},[a(1),a(2),a(3),2]);
B.cores{d,1} = reshape(B.cores{d,1},[a(1),a(2),a(3)*2,a(4)]);
B.dim(d,:) = [a(1),a(2),a(3)*2,a(4)];
% TTm -> TT
[BTT,nrB,ncB] = TTm2TT(B);

%% [A 0] + [0 B]
% check whether the dimensions are the same
if ~isequal(nrA,nrB) || ~isequal(ncA,ncB)
   error('The dimensions the two TT are not identical') 
end
% addTT is done with TT and not with TTm
%% Change format from here on!
CTT = addTT(ATT,BTT);
C = TT2TTm(CTT,nrA,ncA);

%% checked 01/05/2021
end