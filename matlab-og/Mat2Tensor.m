function Tensor = Mat2Tensor(Mat,Rdim,Cdim)
% This function rewrites any matrix NxM into a tensor N1M1xN2M2x...xNnMn, which
% dimensions can be specified by dim.
% INPUT: Matrix (NxM), Row partitioning (N1,N2,...,Nn), Column partitioning (M1,M2,...,Mn)
% OUTPUT: Tensor (N1M1xN2M2x...xNnMn)

%% start function
[Rmat,Cmat] = size(Mat);
if Rmat ~= prod(Rdim)
    error('The chosen dimensions for the rows are not suitable. Please choose differently.')
end
if Cmat ~= prod(Cdim)
    error('The chosen dimensions for the columns are not suitable. Please choose differently.')
end
r = length(Rdim);
c = length(Cdim);
Tensor = reshape(Mat,[Rdim,Cdim]);
order(1:2:r+c) = 1:r;
order(2:2:r+c) = 1+r:c+r;
Tensor = permute(Tensor,order);
Tensor = reshape(Tensor,Rdim.*Cdim);
%% checked on 30/10/2020 & 01/05/2021
end