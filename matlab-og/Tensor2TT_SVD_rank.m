function [TT] = Tensor2TT_SVD_rank(Tensor,r)
% This function computes the Tensor-Train cores and the corresponding ranks
% of a given tensor and epsilon. The SVD approach is taken to compute the
% cores of the TT. The cores are truncated via the rank. The rank can
% therefore be specified beforehand to obtain the desired core dimensions.

%% start function
% d equals the amount of cores
d = length(size(Tensor));
% set C equal to the tensor
C = Tensor; 
% the dimension of the tensor
n = size(Tensor);

% performing SVD from left to right
for k = 2:d 
    % reshape C every time by taking out its left rank and dimension n
    C = reshape(C,[r(k-1)*n(k-1), numel(C)/r(k-1)/n(k-1)]);
    % taking the SVD of C
    [U,S,V] = svd(C,'econ');
    % truncation of U,S,V via the rank
    S_trun = S(1:r(k),1:r(k));
    U_trun = U(:,1:r(k));
    V_trun = V(:,1:r(k));
    % compute G via the truncated U, according to ranks and dim
    G{k-1,1} = reshape(U_trun,[r(k-1),n(k-1),r(k)]);
    % set C equal to truncated SV', as preparation for new iteration
    C = S_trun*V_trun';
end
% final core is equal to the leftover truncated SV'
G{d,1} = C;
TT = G;

%% checked via TT2Tensor.m on 30/10/2020
% note: error is relatively large for some instances
end