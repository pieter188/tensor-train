%% Modified Gram-Schmidt algorithm for 20 by 20 matrix
% Pieter van Klaveren
% 01-04-2021
clear all; clc; close all

A = rand(40,40)*5;
Row = [2,2,5,2];
Col = [2,2,5,2];
Aten = Mat2Tensor(A,Row,Col);
eps = 1e-6;
ATT = Tensor2TT_SVD_eps(Aten,eps);
ATTm = TT2TTm(ATT,Row,Col);
clear Aten ATT
for i = 1:prod(Col)
    % take ith vector out of TTm of A
    % function to determine the code of indexing.
    sizeATTm = size(ATTm);
    e = multi_index(Col,i);
    % column pick from ATTm via e1,e2,e3
    [d,~] = size(ATTm);
    for k = 1:d
        a = size(ATTm{k,1});
        if k == d
           a = [a,1]; 
        end
        V{1,i}{k,1} = reshape(permute(ATTm{k,1},[1,2,4,3]),[numel(ATTm{k,1})/a(3),a(3)])*e{k,1};
        V{1,i}{k,1} = reshape(V{1,i}{k,1},[a(1),a(2),a(4)]);
    end
end   

%% combine the vectors in V again
Vf = combineTT(V,Col);
[Vf,R,C] = TTm2TT(Vf);

%% check if Vf is somewhat equal to A
Vf = TT2Tensor(Vf);
Vf = Tensor2Mat(Vf,R,C);