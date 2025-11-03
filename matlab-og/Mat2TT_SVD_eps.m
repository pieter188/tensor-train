function [TT,r] = Mat2TT_SVD_eps(mat,dim,eps)
% This function computes the Tensor-Train cores and the corresponding ranks
% of a given matrix, desired dimension and epsilon. The SVD approach is taken to compute the
% cores of the TT. The cores are truncated via epsilon, this variable can
% be adjusteda and has direct influence on the truncation.

%% start function
% d equals the amount of cores
d = length(dim);
% delta determines where to truncate in the S matrix
delta = eps/sqrt(d-1)*sqrt(sumsqr(mat));
% set C equal to the tensor
C = mat; 
% r is a vector containing the ranks, to be filled in during iteration
r = zeros(1,d+1);r(1) = 1;r(end) = 1;
% the dimension of the tensor
n = dim;

% performing SVD from left to right
for k = 2:d
    % reshape C every time by taking out its left rank and dimension n
    C = reshape(C,[r(k-1)*n(k-1), numel(C)/r(k-1)/n(k-1)]);
    % taking the SVD of C
    [U,S,V] = svd(C,'econ');
    % delta-truncate SVD, take frobinius norm of singular values, E is what
    % you throw away.
    for i = min(size(S)):-1:1
        a = sqrt(sumsqr(S(i:end,i:end)));
        if a>delta
            i=i;
            break;
        end
    end
    % truncate U,S,V accordingly
    S_trun = S(1:i,1:i);        
    U_trun = U(:,1:i);          
    V_trun = V(:,1:i);          
    % set C equal to truncated USV'
    C = U_trun*S_trun*V_trun';
    % compute rank of C, this will be the right rank of the core
    r(k) = rank(C);
    % compute G via the truncated U, according to ranks and dim
    G{k-1,1} = reshape(U_trun,[r(k-1),n(k-1),r(k)]);
    % set C equal to truncated SV', as preparation for new iteration
    C = S_trun*V_trun';
end
% final core is equal to the leftover truncated SV'
G{d,1} = C;
TT = G;

%% checked via reconstruction on 31/10/2020 AND via Tensor2TT_SVD_eps.m on 31/10/2020
end
