function [TT,Rdim,Cdim] = TTm2TT(TTm)
% This function transforms a TensorTrain matrix into a TensorTrain. It will
% do so by combining the dimensions of the TTm-cores into one TT-core. Thus
% the cores will go from being a 4-way tensor to a 3-way tensor.
% INPUT: TTm (Tensor Train matrix)
% OUTPUT: TT (Tensor Train), Rdim (Row dimensions), Cdim (Column dimensions)

% d equals the amount of cores of TTm
d = size(TTm.dim,1);
for i = 1:d
    % obtain TT-core by reshaping TTm-core
    TT.cores{i,1} = reshape(TTm.cores{i,1},[TTm.dim(i,1),TTm.dim(i,2)...
        *TTm.dim(i,3),TTm.dim(i,4)]);
    TT.dim(i,:) = [TTm.dim(i,1),TTm.dim(i,2)*TTm.dim(i,3),TTm.dim(i,4)];
    Rdim(i) = TTm.dim(i,2);
    Cdim(i) = TTm.dim(i,3);
end
%% checked on 01/05/2021
end

% %% start function
% % d equals the amount of cores of TTm
% [d,~] = size(TTm);
% for i = 1:d
%     if i == 1
%         % obtain size of each TTm-core
%         [a1,Rdim(i),Cdim(i),a4] = [size(TTm{i,1})];
%         % obtain TT-core by reshaping TTm-core
%         TT{i,1} = reshape(TTm{i,1},[a1,Rdim(i)*Cdim(i),a4]);
%     elseif i == d
%         % obtain size of each TTm-core
%         [a1,Rdim(i),Cdim(i)] = [size(TTm{i,1})];
%         % obtain TT-core by reshaping TTm-core
%         TT{i,1} = reshape(TTm{i,1},[a1,Rdim(i)*Cdim(i)]);
%     else
%         % obtain size of each TTm-core
%         [a1,Rdim(i),Cdim(i),a4] = size(TTm{i,1});
%         % obtain TT-core by reshaping TTm-core
%         TT{i,1} = reshape(TTm{i,1},[a1,Rdim(i)*Cdim(i),a4]);
%     end
% end
% %% checked on 30/10/2020