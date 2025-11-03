function mat = Tensor2Mat_moden(tensor,mode)
y = size(tensor);
d = length(y);
% x = ceil(input(['What mode-n matricization do you want? please provide a number between 1 and ' num2str(d) '? ']));
if mode > d
    disp('This mode-n matricization is not possible, please try again')
    mode = ceil(input(['What mode-n matricization do you want? please provide a number between 1 and ' num2str(d) '? ']));
end
    

while d>=4 % in this if-statement I set the mode-n at the front of the matrix
    switch mode
        case 1
            mat = reshape(permute(tensor,[mode,mode+1:d]),[y(mode),y(mode+1:d-2),y(d-1)*y(d),1]);
        case 2
            mat = reshape(permute(tensor,[mode,1,mode+1:d]),[y(mode),y(1),y(mode+1:d-2),y(d-1)*y(d),1]);
        case d-1
            mat = reshape(permute(tensor,[mode,1:mode-1,d]),[y(mode),y(1:mode-2),y(mode-1)*y(d),1]);
        case d
            mat = reshape(permute(tensor,[mode,1:mode-1]),[y(mode),y(1:mode-3),y(mode-2)*y(mode-1),1]);
        otherwise
            mat = reshape(permute(tensor,[mode,1:mode-1,mode+1:d]),[y(mode),y(1:mode-1),y(mode+1:d-2),y(d-1)*y(d),1]);
    end
    y = size(mat); % update length y
    d = length(y); % substract 1 from d
    mode = 1; % since right fiber is at the front of the matrix already we can set x to 1
    tensor = mat;
end 

if d == 3
    switch mode
        case 1
            mat = reshape(permute(tensor,[mode,mode+1:d]),[y(mode),y(d-1)*y(d),1]);
        case 2
            mat = reshape(permute(tensor,[mode,1,d]),[y(mode),y(1)*y(d),1]);
        case 3
            mat = reshape(permute(tensor,[mode,1:mode-1]),[y(mode),y(1)*y(2),1]);
    end
end

end













%% for-loop style
% case 1
%         disp('mode-1 selected') % mode 1: column-fibers of tensor are columns of matrix
%         for k = 1:z
%             for j = 1:y
%                A_mat = [A_mat,A(:,j,k)];
%             end
%         end
%     case 2
%         disp('mode-2 selected') % mode 2: row-fibers of tensor are columns of matrix
%         for k = 1:z
%             for i = 1:x
%                A_mat = [A_mat,A(i,:,k)'];
%             end
%         end
%     case 3
%         disp('mode-3 selected') % mode 3: tube-fibers of tensor are columns of matrix
%         for j = 1:y
%             for i = 1:x
%                A_mat = [A_mat,permute(A(i,j,:),[3,1,2])];
%             end
%         end