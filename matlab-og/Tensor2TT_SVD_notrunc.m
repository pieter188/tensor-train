function [TT,r] = Tensor2TT_SVD_notrunc(Tensor)
% This function computes the Tensor-Train cores and the corresponding ranks
% of a given tensor and epsilon. The SVD approach is taken to compute the
% cores of the TT. The cores are not truncated.  

%% start function
% d equals the amount of cores
d = length(size(Tensor));
% set C equal to the tensor
C = Tensor; 
% r is a vector containing the ranks, to be filled in during iteration
r = zeros(1,d+1);r(1) = 1;r(end) = 1;
% the dimension of the tensor
n = size(Tensor);

% performing SVD from left to right
for k = 2:d
    % reshape C every time by taking out its left rank and dimension n
    C = reshape(C,[r(k-1)*n(k-1), numel(C)/r(k-1)/n(k-1)]);
    % taking the SVD of C
    [U,S,V] = svd(C,'econ');
    % set C equal to USV'
    C = U*S*V';
    % compute rank of C, this will be the right rank of the core
    r(k) = rank(C);
    % compute G via the truncated U, according to ranks and dim
    G{k-1,1} = reshape(U,[r(k-1),n(k-1),r(k)]);
    % set C equal to truncated SV', as preparation for new iteration
    C = S*V';
end
% final core is equal to the leftover truncated SV'
G{d,1} = C;
TT = G;

%% checked via TT2Tensor.m on 30/10/2020
end
