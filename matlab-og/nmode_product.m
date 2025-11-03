function A_n = nmode_product(Tensor,mat,mode)
% this function computes the n-mode product of a matrix/vector times a
% tensor. This is done by using the matricization.m file to obtain the
% unfolded tensor in matrix format, then simply multiplying the
% matrix/vector with the unfolded tensor, and then folding the tensor up
% again (in reverse order w.r.t. the matricization) and lastly using the
% permute function to obtain the right order of the dimensions.

% X * U with dim(X) = I1*...*IN and dim(U) = J*In becomes dim(X*U) =
% i1*..*In-1*J*In+1*..*IN

y = size(Tensor);
d = length(y);
% x = ceil(input(['What n-mode product do you want? please provide a number between 1 and ' num2str(d) '? ']));

% check size of matrix/vector U
[aa, bb]=size(mat);
if bb ~= y(mode)
    disp('matrix dimensions do not match')
end
    
% use matricization from previous question
A_mat = Tensor2Mat_moden(Tensor,mode);

A_n = mat*A_mat; % atm: dim(A_n) = Jx(I1*...*In-1*In+1*...*IN)
% now fold new matrix again to dim(A_n) = I1x...xIn-1xJxIn+1x...xIN
% information about the original matrix a is stored in y and d

y1 = [y(1:mode-1),y(mode+1:d)]; % dimensions of A without dim n, used in forloop for folding the matrix
y_mod = size(A_n);      % gives the curring size of the to-be folded matrix
cs = length(y_mod);     % gives the amount of dimensions of the to be folded matrix
p = 1; %parameter to know over which axis to fold
while p <= d-1
    A_n = reshape(A_n,[y_mod(1:end-1),y1(p),y_mod(end)/y1(p)]);   % folding from 1->d
    y_mod = size(A_n);
    cs = length(y_mod);
    p = p+1;
end
A_n = permute(A_n,[2:mode,1,mode+1:d]);
end






%%
% while d <= 3
%     switch x
%         case 1
%             A_n = permute(reshape(A_n,[aa,y(x+1:end)]),[1,2,3]);
%         case 2
%             A_n = permute(reshape(A_n,[aa,y(1:x-1),y(x+1:end)]),[2,1,3]);
%         case 3
%             A_n = permute(reshape(A_n,[aa,y(1:x-1)]),[2,3,1]);
%     end
% end
