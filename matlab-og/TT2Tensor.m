function Tensor = TT2Tensor(TT)
% This function reconstructs a TT into a tensor by performing contractions 
% over the TT-ranks.
% INPUT: TT
% OUTPUT: Tensor

%% start function
% d is the amount of TT cores
d = size(TT.dim,1);   
% size of the first TT core
a = TT.dim(1,:);  
% save i_1 in separate vector
c(1) = a(2);   
% reshape so that rank is 'free'
Tensor = reshape(TT.cores{1,1},[a(1)*a(2),a(3)]);     
for i = 2:d
    % size of the ith TT core
    a = TT.dim(i,:);
    % save ith dim in separate vector
    c(i) = a(2);     
    % size of B_rec
    b = size(Tensor);                
    % multiply B_rec with the ith core, (summation over rank)
    Tensor = Tensor*reshape(TT.cores{i,1},[a(1),a(2)*a(3)]);   
    % reshape to make the rank 'free' again for the next multiplication
    Tensor = reshape(Tensor,[b(1)*a(2),a(3)]);            
end
% reconstruct tensor using c
Tensor = reshape(Tensor,c);           

%% Checked on 04/05/2021
end

% % d is the amount of TT cores
% [d,~] = size(TT);   
% % size of the first TT core
% a = size(TT{1,1});  
% % save i_1 in separate vector
% c(1) = a(2);   
% % reshape so that rank is 'free'
% Tensor = reshape(TT{1,1},[a(1)*a(2),a(3)]);     
% for i = 2:d
%     % size of the ith TT core
%     a = size(TT{i,1});  
%     % save ith dim in separate vector
%     c(i) = a(2);     
%     % size of B_rec
%     b = size(Tensor);                
%     if i == d
%         a = [size(TT{i,1}),1];
%     end
%     % multiply B_rec with the ith core, (summation over rank)
%     Tensor = Tensor*reshape(TT{i,1},[a(1),a(2)*a(3)]);   
%     % reshape to make the rank 'free' again for the next multiplication
%     Tensor = reshape(Tensor,[b(1)*a(2),a(3)]);            
% end
% % reconstruct tensor using c
% Tensor = reshape(Tensor,c);           
% 
% %% Checked via Tensor2TT_SVD_eps on 30/10/2020 and on 26/03/2021