function Tensor = TTm2Tensor(TTm)
% This function reconstructs TT-format into a tensor.
%% start function
[d,~] = size(TTm);
% convert TTm-core to TT-core, thus the dimensions of the TT are
% Rdim(i)*Cdim(i)
[TT,Rdim,Cdim] = TTm2TT(TTm);
% tensor dimensions are Rdim(i)*Cdim(i)
Tensor = TT2Tensor(TT);
% split row and columns
dim(1:2:2*d-1) = Rdim;
dim(2:2:2*d) = Cdim;
Tensor = reshape(Tensor,dim);
% reshape tensor to obtain row first and then columns
Tensor = permute(Tensor,[1:2:2*d-1,2:2:2*d]);

%% Checked on 31/10/2020 via TTm2TT.m and TT2Tensor.m and some reshaping.
end




% 
% % d is the amount of TT cores
% [d,~] = size(TTm); 
% % retrieve size of TTm-core
% [a1,Rdim(1),Cdim(1),a4] = size(TTm{1,1});
% % reshape so that rank is 'free'
% mat = reshape(TTm{1,1},[a1*Rdim(1)*Cdim(1),a4]);     
% for i = 2:d
%     % retrieve size of TTm-core
%     [a1,Rdim(i),Cdim(i),a4] = size(TTm{i,1});
%     % size of B_rec
%     [Rmat,~] = size(mat);                
%     % multiply B_rec with the ith core, (summation over rank)
%     mat = mat*reshape(TTm{i,1},[a1,Rdim(i)*Cdim(i)*a4]);
%     % reshape to make the rank 'free' again for the next multiplication
%     mat = reshape(mat,[Rmat*Rdim(i)*Cdim(i),a4]);            
% end
% % at the end of the for-loop mat has the row and column dimensions
% % alternating
% order = 
% 
% % reconstruct tensor using c
% Tensor = reshape(mat,order);    
% % permute to have all the upper i's of the cores first and then all the under i's
% Tensor = permute(Tensor,[1:2:2*d-1,2:2:2*d]); 