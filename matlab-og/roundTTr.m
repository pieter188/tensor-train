function [TTt,rt] = roundTTr(TT,r)
% start TT rounding algorithm, this function must truncate the TT-ranks of
% the tensor Train of A such that the relative error is less than epsilon.
% This is realized by taking the RQ-factorization from one side to the
% other side and then taking the SVD-decomposition and truncating by preset
% ranks.

d = size(TT.dim,1);

TTn = sitek(TT,d);
norm2 = norm(reshape(TTn.cores{d,1},[numel(TTn.cores{d,1}),1]),'fro');
delta = eps/sqrt(d-1)*norm2;


% the norm of the tensor train is now in the last tensor train matrix.
%% checked 04/05/2021
% right to left SVD-decomposition procedure
n = TTn.dim(:,2)';
[U,S,V] = svd(TTn.cores{d,1},'econ');

S_trun = S(1:r(d),1:r(d));
U_trun = U(:,1:r(d));
V_trun = V(:,1:r(d));

TTt.cores{d,1} = reshape(V_trun',[r(d),n(d),r(d+1)]);
TTt.dim(d,:) = [r(d),n(d),r(d+1)];
C = U_trun*S_trun;
for j = d-1:-1:2
    a = TTn.dim(j,:);
    b = size(C);
    TTn.cores{j,1} = reshape(reshape(TTn.cores{j,1},[a(1)*a(2),a(3)])*C,...
        [a(1),a(2)*b(2)]);
    [U,S,V] = svd(TTn.cores{j,1},'econ');
    S_trun = S(1:r(j),1:r(j));
    U_trun = U(:,1:r(j));
    V_trun = V(:,1:r(j));
    TTt.cores{j,1} = reshape(V_trun',[r(j),n(j),r(j+1)]);
    TTt.dim(j,:) = [r(j),n(j),r(j+1)];
    C = U_trun*S_trun;
end
a = TTn.dim(1,:);
TTt.cores{1,1} = permute(reshape(TTn.cores{1,1},[a(1)*a(2),a(3)])*C,[3,1,2]);
TTt.dim(1,:) = [1,a(1)*a(2),size(C,2)];

end


% [d,~] = size(TT);
% % left to right QR-decomposition procedure to bring the norm to the right
% % most core
% a = size(TT{1,1});
% % combine left-rank and index
% TT{1,1} = reshape(TT{1,1},[a(1)*a(2),a(3)]);
% % qr-decomposition of first core
% [Q,R] = qr(TT{1,1});
% % reshape Q, orthonormal matrix, into new tensor train core
% TTn{1,1} = reshape(Q,[a(1),a(2),a(1)*a(2)]);
% for i = 2:d-1
%     a = size(TT{i,1});
%     % obtain size of R
%     b = size(R);
%     % combine the norm of the previous core, R with the current core TT{i,1}
%     TT{i,1} = R*reshape(TT{i,1},[a(1),a(2)*a(3)]);
%     TT{i,1} = reshape(TT{i,1},[b(1),a(2),a(3)]);
%     TT{i,1} = reshape(TT{i,1},[b(1)*a(2),a(3)]);
%     [Q,R] = qr(TT{i,1});
%     TTn{i,1} = reshape(Q,[b(1),a(2),b(1)*a(2)]);
% end
% a = size(TT{d,1});
% % obtain size of R
% b = size(R);
% % combine the norm of the previous core, R with the current core TT{i,1}
% TTn{d,1} = R*TT{d,1};
% % the norm of the tensor train is now in the last tensor train matrix.
% 
% % right to left SVD-decomposition procedure
% n = zeros(d,1);
% a = size(TTn{d,1});
% n(d) = a(2);
% [U,S,V] = svd(TTn{d,1},'econ');
% S_trun = S(1:r(d),1:r(d));
% U_trun = U(:,1:r(d));
% V_trun = V(:,1:r(d));
% 
% C = U_trun*S_trun*V_trun';
% rt(d) = rank(C);
% % if rt(d) ~= r(d)
% %     error;
% % end
% TTt{d,1} = reshape(V_trun',[r(d),n(d),r(d+1)]); 
% C = U_trun*S_trun;
% for j = d-1:-1:2
%     a = size(TTn{j,1});
%     n(j) = a(2);
%     b = size(C);
%     [U,S,V] = svd(reshape(reshape(TTn{j,1},[a(1)*a(2),a(3)])*C,[a(1),a(2)*b(2)]),'econ');
%     S_trun = S(1:r(d),1:r(d));
%     U_trun = U(:,1:r(d));
%     V_trun = V(:,1:r(d));
%     C = U_trun*S_trun*V_trun';
%     rt(j) = rank(C);
% %     if rt(j) ~= r(j)
% %         error;
% %     end
%     TTt{j,1} = reshape(V_trun',[r(j),n(j),r(j+1)]); 
%     C = U_trun*S_trun;
% end
% a = size(TTn{1,1});
% TTt{1,1} = permute(reshape(TTn{1,1},[a(1)*a(2),a(3)])*C,[3,1,2]);