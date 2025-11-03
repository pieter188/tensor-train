function G = innerproduct(A,B)
% computes the inner products of two tensors

%check size of tensors
[x] = size(A);
[y] = size(B);
if x==y
    G = times(A,B);
else
    fprintf('tensors do not have the same size')
end
G = sumsqr(G);
end