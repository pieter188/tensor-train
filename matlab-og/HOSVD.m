function [G,A] = HOSVD(I,R)
% This function computes Tuckers first method (Higher-Order SVD).
    d=length(size(I));
    
    for i = 1:d
        X = Tensor2Mat_moden(I,i);
        [U,S,V] = svd(X,'econ');
        if R(i)>max(size(U))
            R(i) = max(size(U));
        end
        A{1,i} = U(:,1:R(i));
    end
    G = I;
    for j = 1:d
        G = nmode_product(G,A{1,j}',j);
    end
    HOSVD_storage = numel(G)+numel(A{1,1})+numel(A{1,2})+numel(A{1,3})+numel(A{1,4})+numel(A{1,5})+numel(A{1,6})+numel(A{1,7})+numel(A{1,8})+numel(A{1,9});
end