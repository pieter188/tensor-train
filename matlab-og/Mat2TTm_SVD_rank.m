function TTM = Mat2TTm_SVD_rank(mat,dim,r)
% This function computes a Tensor Train of a given matrix A, with a given
% desired dimension of the TT-cores dim and a relative approximation error.
if numel(mat) ~= prod(dim)
    error('The chosen dimensions are not suitable. Please choose differently.')
end
d = length(dim)/2;
C = mat; 
n = dim;
% only useful to use svd for low rank.

for k = 2:d % here I compute B = G1(r1,i1,r2)*G(r2,i2,r3)*...
    C = reshape(C,[r(k-1)*n(k-1)*n(d+k-1), numel(C)/r(k-1)/n(k-1)/n(d+k-1)]);
    [U,S,V] = svd(C,'econ');

    S_trun = S(1:r(k),1:r(k));
    U_trun = U(:,1:r(k));
    V_trun = V(:,1:r(k));
    
    G{k-1,1} = reshape(U_trun,[r(k-1),n(k-1),n(d+k-1),r(k)]);
    C = S_trun*V_trun';
end
G{d,1} = reshape(C,[r(d),n(d),n(d+d),r(d+1)]);
TTM = G;
end