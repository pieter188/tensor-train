function C = concatTT(A,B)
% input: A and B, two TTm's
% output: C, TT concatenation of A and B
% In the first section this function creates two structures [A 0] and [0
% B]. It is required that A and B have the same dimension in order to
% succesfully do the concatenation.
if size(A) ~= size(B)
    error('The matrix dimensions of A and B are not the same')
end 

%% [A 0]
[d,~] = size(A);
V = [1 0];
% connect V to the last core of A
a = [size(A{d,1}),1];
A{d,1} = reshape(A{d,1},[a(1)*a(2),a(3)]);
A{d,1} = A{d,1}*V;
A{d,1} = reshape(A{d,1},[a(1),a(2),2]);
A{d,1} = reshape(A{d,1},[a(1),a(2),2,a(3)]);

% TTm -> TT
[ATT,nrA,ncA] = TTm2TT(A);

%% [0 B]
W = [0 1];
% connect V to the last core of A
a = [size(B{d,1}),1];
B{d,1} = reshape(B{d,1},[a(1)*a(2),a(3)]);
B{d,1} = B{d,1}*W;
B{d,1} = reshape(B{d,1},[a(1),a(2),2]);
B{d,1} = reshape(B{d,1},[a(1),a(2),2,a(3)]);

% TTm -> TT
[BTT,nrB,ncB] = TTm2TT(B);

%% [A 0] + [0 B]
% check whether the dimensions are the same
if ~isequal(nrA,nrB) || ~isequal(ncA,ncB)
   error('The dimensions the two TT are not identical') 
end
% addTT is done with TT and not with TTm
CTT = addTT(ATT,BTT);
CTT = sitek(CTT,1);
C = TT2TTm(CTT,nrA,ncA);

end