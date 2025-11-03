function TTc = sitek(TT,k)
% This function brings any TensorTrain into site-k mixed canonical form.
% The norm of this reshaped tensor train is all in the k-th TT-core.
% INPUT: Tensor-Train, value between 1 and the total amount of cores
% OUTPUT: k-Orthogonal Tensor Train

d = size(TT.dim,1);
if k == 1
    % collect the norm of the dth to the k+1th core
    a = TT.dim(d,:);
    [U,S,V] = svd(TT.cores{d,1},'econ');
    TTc.cores{d,1} = V';
    TTc.dim(d,:) = [size(V'),1];
    aSo = size(S);
    for i = d-1:-1:k+1
        a = TT.dim(i,:);
        [U,S,V] = svd(reshape(reshape(TT.cores{i,1},[a(1)*a(2),a(3)])*U*S,[a(1),a(2)*aSo(2)]),'econ');
        aSn = size(S);
        TTc.cores{i,1} = reshape(V',[aSn(1),a(2),aSo(2)]);
        TTc.dim(i,:) = [aSn(1),a(2),aSo(2)];
        aSo = aSn;
        if i == k+1
           C_right = U*S; 
        end
    end
    a = TT.dim(k,:);
    aSn = size(C_right);
    TTc.cores{k,1} = reshape(reshape(TT.cores{k,1},[a(1)*a(2),a(3)])*C_right,[a(1),a(2),aSn(2)]);
    TTc.dim(k,:) = [a(1),a(2),aSn(2)];
end

if k == d
    % collect the norm of the 1st to the d-1th core
    a = TT.dim(1,:);
    [U,S,V] = svd(reshape(TT.cores{1,1},[a(1)*a(2),a(3)]),'econ');
    aSo = size(S);
    TTc.cores{1,1} = reshape(U,[a(1),a(2),aSo(1)]);
    TTc.dim(1,:) = [a(1),a(2),aSo(1)];
    for i = 2:1:d-1
        a = TT.dim(i,:);
        [U,S,V] = svd(reshape(S*V'*reshape(TT.cores{i,1},[a(1),a(2)*a(3)]),[aSo(1)*a(2),a(3)]),'econ');
        aSn = size(S);
        TTc.cores{i,1} = reshape(U,[aSo(1),a(2),aSn(2)]);
        TTc.dim(i,:) = [aSo(1),a(2),aSn(2)];
        aSo = aSn;
        if i == d-1
           C_left = S*V'; 
        end
    end
    a = TT.dim(k,:);
    aSn = size(C_left);
    TTc.cores{k,1} = reshape(C_left*reshape(TT.cores{k,1},[a(1),a(2)*a(3)]),[aSn(1),a(2),a(3)]);
    TTc.dim(k,:) = [aSn(1),a(2),a(3)];
end
%% checked for case k==1 & k==d on 05/05/2021
%% collect the norm from the 1st to the k-1th core

% a = size(TT.dim,1);
% [U,S,V] = svd(reshape(permute(TT{1,1},[2,3,1]),[a(1)*a(2),a(3)]),'econ');
% TTc{1,1} = permute(U,[3,1,2]);
% for i = 2:k-1
%     a = size(TT{i,1});
%     [U,S,V] = svd(reshape(S*V'*reshape(permute(TT{i,1},[1,2,3]),[a(1),a(2)*a(3)]),[a(1)*a(2),a(3)]),'econ');
%     TTc{i,1} = reshape(U,[a(1),a(2),a(3)]);
%     if i ==k-1
%        C_left = S*V'; 
%     end   
% end

%% collect the norm of the dth to the k+1th core
% a = TT.dim(d,:);
% [U,S,V] = svd(TT.cores{d,1},'econ');
% TTc.cores{d,1} = V';
% TTc.dim(d,:) = [size(V'),1];
% aSo = size(S);
% for i = d-1:-1:k+1
%     a = TT.dim(i,:);
%     [U,S,V] = svd(reshape(reshape(TT.cores{i,1},[a(1)*a(2),a(3)])*U*S,[a(1),a(2)*aSo(2)]),'econ');
%     aSn = size(S);
%     TTc.cores{i,1} = reshape(V',[aSn(1),a(2),aSo(2)]);
%     TTc.dim(i,:) = [aSn(1),a(2),aSo(2)];
%     aSo = aSn;
%     if i == k+1
%        C_right = U*S; 
%     end   
% end
% 
% %% ATT_new(k) equals C_left*ATT(k)*C_right
% % if k ~= 1 && k ~= d
% %     a = size(TT{k,1});
% %     TTc{k,1} = reshape(reshape(reshape(C_left*reshape(TT{k},[a(1),a(2)*a(3)]),[a(1),a(2),a(3)]),[a(1)*a(2),a(3)])*C_right,[a(1),a(2),a(3)]);
% % end
% if k == 1
%     a = TT.dim(k,:);
%     aSn = size(C_right);
%     TTc.cores{k,1} = reshape(reshape(TT.cores{k,1},[a(1)*a(2),a(3)])*C_right,[a(1),a(2),aSn(2)]);
%     TTc.dim(k,:) = [a(1),a(2),aSn(2)];
% end
%% checked 03/05/2021
end

% % collect the norm from the 1st to the k-1th core
% [d,~] = size(TT);
% % a = size(TT{1,1});
% % [U,S,V] = svd(reshape(permute(TT{1,1},[2,3,1]),[a(1)*a(2),a(3)]),'econ');
% % TTc{1,1} = permute(U,[3,1,2]);
% % for i = 2:k-1
% %     a = size(TT{i,1});
% %     [U,S,V] = svd(reshape(S*V'*reshape(permute(TT{i,1},[1,2,3]),[a(1),a(2)*a(3)]),[a(1)*a(2),a(3)]),'econ');
% %     TTc{i,1} = reshape(U,[a(1),a(2),a(3)]);
% %     if i ==k-1
% %        C_left = S*V'; 
% %     end   
% % end
% 
% %% collect the norm of the dth to the k+1th core
% a = size(TT{d,1});
% [U,S,V] = svd(TT{d,1},'econ');
% TTc{d,1} = V';
% aSo = size(S);
% for i = d-1:-1:k+1
%     a = size(TT{i,1});
%     if length(a) == 2
%         a = [a,1];
%     end
%     [U,S,V] = svd(reshape(reshape(TT{i,1},[a(1)*a(2),a(3)])*U*S,[a(1),a(2)*aSo(2)]),'econ');
%     aSn = size(S);
%     TTc{i,1} = reshape(V',[aSn(1),a(2),aSo(2)]);
%     aSo = aSn;
%     if i == k+1
%        C_right = U*S; 
%     end   
% end
% 
% %% ATT_new(k) equals C_left*ATT(k)*C_right
% if k ~= 1 && k ~= d
%     a = size(TT{k,1});
%     TTc{k,1} = reshape(reshape(reshape(C_left*reshape(TT{k},[a(1),a(2)*a(3)]),[a(1),a(2),a(3)]),[a(1)*a(2),a(3)])*C_right,[a(1),a(2),a(3)]);
% end
% if k == 1
%     a = size(TT{k,1});
%     if length(a) == 2
%         a = [a,1];
%     end
%     aSn = size(C_right);
%     TTc{k,1} = reshape(reshape(TT{k,1},[a(1)*a(2),a(3)])*C_right,[a(1),a(2),aSn(2)]);
% end