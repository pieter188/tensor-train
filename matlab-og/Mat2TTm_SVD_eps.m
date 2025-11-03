function TTm = Mat2TTm_SVD_eps(mat,Rdim,Cdim,eps)
% This function computes a Tensor Train of a given matrix A, with a given
% desired dimension of the TT-cores dim and a relative approximation error.

%% start function
% create tensor out of matrx with specified dimensions for row and column
Tensor = Mat2Tensor(mat,Rdim,Cdim);
% amount of cores in TT(m)
d = length(Rdim);
% permute matrix to obtain right order of alternating row and column
% dimensions
orderODD = 1:2:2*d-1;
orderEVEN = 2:2:2*d;
order(orderODD) = 1:d;
order(orderEVEN) = d+1:2*d;
Tensor = permute(Tensor,order);
% constract each corresponding row and column dimension
Tensor = reshape(Tensor,[Rdim.*Cdim]);

% create TensorTrain from Tensor
TT = Tensor2TT_SVD_eps(Tensor,eps);

% reshape TT into TTm
TTm = TT2TTm(TT,Rdim,Cdim);

%% checked on 31/10/2020
end