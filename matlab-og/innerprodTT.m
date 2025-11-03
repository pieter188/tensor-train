function [c] = innerprodTT(TT1,TT2)
% This function computes the innerproduct of two TT
d = size(TT1.dim,1);

for i = 1:d
    a = TT1.dim(i,:);
    b = TT2.dim(i,:);
%     TT.cores{i,1} = reshape(permute(TT1.cores{i,1},[1,3,2]),...
%         [numel(TT1.cores{i,1})/a(2),a(2)])*reshape(permute(TT2.cores{i,1}...
%         ,[2,1,3]),[b(2),numel(TT2.cores{i,1})/b(2)]);
%     if i ~= 1 || i ~= d
    TT.cores{i,1} = reshape(permute(reshape(reshape(permute(TT1.cores{i,1},...
        [1,3,2]),[numel(TT1.cores{i,1})/a(2),a(2)])*reshape(permute...
        (TT2.cores{i,1},[2,1,3]),[b(2),numel(TT2.cores{i,1})/b(2)]),...
        [a(1),a(3),b(1),b(3)]),[1,3,2,4]),[a(1)*b(1),a(3)*b(3)]);
%     end
end
c = TT.cores{1,1};
for j = 2:d
    c = c*TT.cores{j,1};    
end
%% checked 03/05/2021
end