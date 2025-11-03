function TTm = TT2TTm(TT,Rdim,Cdim)
% This function transforms a TensorTrain into a TensorTrain matrix. It will
% do so by splitting the dimension of the TT-cores into two dimensions of 
% TT-cores. This will be done so that the row dimensions of the original
% matrix comes before the column dimension of the original matrix.
% INPUT: TT (Tensor Train), Rdim (Row dimensions), Cdim (Column dimensions)
% OUTPUT: TTm (Tensor Train matrix)


%% start function
% d equals the amount of cores of TTm
d = size(TT.dim,1);
for i = 1:d
    % obtain TTm-core by reshaping the TT-core
    TTm.cores{i,1} = reshape(TT.cores{i,1},[TT.dim(i,1),Rdim(i),Cdim(i),TT.dim(i,3)]);
    % obtain TTm core dimensions
    TTm.dim(i,:) = [TT.dim(i,1),Rdim(i),Cdim(i),TT.dim(i,3)];
end
%% checked on 01/05/2021
end