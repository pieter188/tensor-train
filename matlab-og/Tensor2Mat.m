function Mat = Tensor2Mat(Tensor,Rdim,Cdim)
% This function outputs the original matrix from a given tensor.
% Tensor: R1C1xR2C2x...xRNCN, Marix: R1*R2*...*RN;C1*C2*...*CN.
% Rdim and Cdim denote the partitioning of row and columns dimensions
% respectively.
% INPUT: Tensor, Rdim (Row dimensions), Cdim (Column dimensions)
% OUTPUT: Matrix

%% start function
% obtain tensor size
dim = size(Tensor);
% check of the desired number of row and column dimensions add up to the total
% amount of dimensions of the tensor
if prod(dim) ~= prod(Rdim)*prod(Cdim)
    error('The number of rows and columns selected for the tensor are not suitable. Please choose differently.')
end
dimOrder = zeros(1,2*length(dim));
dimOrder(1:2:length(dimOrder)) = Rdim;
dimOrder(2:2:length(dimOrder)) = Cdim;
% reshape the tensor according to the row and column dimensions
Mat = reshape(Tensor,dimOrder);
% permute the matrix such that the rows and columns are split:
% R1xR2...xRNxC1xC2...xCN.
Mat = permute(Mat,[1:2:length(dimOrder), 2:2:length(dimOrder)]);
% put the row and column partitionings back together
Mat = reshape(Mat,[prod(Rdim), prod(Cdim)]);

%% checked on 16/03/2021 and 26/03/2021
end