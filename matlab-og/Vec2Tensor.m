function Tensor = Vec2Tensor(Vec,Rdim)
% this function rewrites any matrix NxM into a tensor N1M1xN2M2x...xNnMn, which
% dimensions can be specified by dim.
%% start function
Rmat = length(Vec);
if Rmat ~= prod(Rdim)
    error('The chosen dimensions for the rows are not suitable. Please choose differently.')
end
Tensor = reshape(Vec,Rdim);
%% checked on 30/10/2020 & 01/05/2021
end