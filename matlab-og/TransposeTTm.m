function Transpose = TransposeTTm(TTm)
% This function computes the transpose of a matrix by switching the 2nd and
% 3th index of each of the TTm-cores.
% input: TTm
% output: TTm'

% obtain the amount of cores d.
d = size(TTm.dim,1);
for i = 1:d
    Transpose.cores{i,1} = permute(TTm.cores{i,1},[1,3,2,4]);
    Transpose.dim(i,:) = [TTm.dim(i,1),TTm.dim(i,3),TTm.dim(i,2),TTm.dim(i,4)];
end
%% checked on: 01/05/2021
end


% [d,~] = size(TTm);             % obtain the amount of cores d. 
% for i = 1:d
%     Transpose{i,1} = permute(TTm{i,1},[1,3,2,4]);
% end