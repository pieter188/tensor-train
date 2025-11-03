function [TTt,r] = roundTT(TT,eps)
% start TT rounding algorithm, this function must truncate the TT-ranks of
% the tensor Train of A such that the relative error is less than epsilon.
% This is realized by taking the RQ-factorization from one side to the
% other side and then taking the SVD-decomposition with a delta trunctation
% to obtain the truncated TT of a existing TT.
TTori = TT;
d = size(TT.dim,1);
TTori = sitek(TTori,1);
norm2 = norm(reshape(TTori.cores{1,1},[numel(TTori.cores{1,1}),1]),'fro');
delta = eps/sqrt(d-1)*norm2;
% left to right QR-decomposition procedure to bring the norm to the right
% most core
a = TT.dim(1,:);
% combine left-rank and index
TT.cores{1,1} = reshape(TT.cores{1,1},[a(1)*a(2),a(3)]);
% qr-decomposition of first core
[Q,R] = qr(TT.cores{1,1});
% reshape Q, orthonormal matrix, into new tensor train core
TTn.cores{1,1} = reshape(Q,[a(1),a(2),a(1)*a(2)]);
TTn.dim(1,:) = [a(1),a(2),a(1)*a(2)];
for i = 2:d-1
    a = TT.dim(i,:);
    % obtain size of R
    b = size(R);
    % combine the norm of the previous core, R with the current core TT{i,1}
    TT.cores{i,1} = R*reshape(TT.cores{i,1},[a(1),a(2)*a(3)]);
    TT.cores{i,1} = reshape(TT.cores{i,1},[b(1),a(2),a(3)]);
    TT.cores{i,1} = reshape(TT.cores{i,1},[b(1)*a(2),a(3)]);
    [Q,R] = qr(TT.cores{i,1});
    TTn.cores{i,1} = reshape(Q,[b(1),a(2),b(1)*a(2)]);
    TTn.dim(i,:) = [b(1),a(2),b(1)*a(2)];
end
a = TT.dim(d,:);
% obtain size of R
b = size(R);
% combine the norm of the previous core, R with the current core TT{i,1}
TTn.cores{d,1} = R*TT.cores{d,1};
TTn.dim(d,:) = [b(1),a(2),a(3)];

% the norm of the tensor train is now in the last tensor train matrix.
%% checked 04/05/2021
% right to left SVD-decomposition procedure
r = zeros(1,d+1);r(1) = 1;r(d+1) = 1;
n = TTn.dim(:,2)';
[U,S,V] = svd(TTn.cores{d,1},'econ');

% delta-truncate SVD, take frobinius norm of singular values, E is what
% you throw away.
for k = min(size(S)):-1:1
    a = norm(S(k:min(size(S)),k:min(size(S))),'fro');
    if a>delta
        k = k;
        break;
    end
end
S_trun = S(1:k,1:k);
U_trun = U(:,1:k);
V_trun = V(:,1:k);

C = U_trun*S_trun*V_trun';
r(d) = rank(C);
TTt.cores{d,1} = reshape(V_trun',[r(d),n(d),r(d+1)]);
TTt.dim(d,:) = [r(d),n(d),r(d+1)];
C = U_trun*S_trun;
for j = d-1:-1:2
    a = TTn.dim(j,:);
    b = size(C);
    TTn.cores{j,1} = reshape(reshape(TTn.cores{j,1},[a(1)*a(2),a(3)])*C,...
        [a(1),a(2)*b(2)]);
    [U,S,V] = svd(TTn.cores{j,1},'econ');

    % delta-truncate SVD, take frobinius norm of singular values, E is what
    % you throw away.
    for k = min(size(S)):-1:1
        a = norm(S(k:min(size(S)),k:min(size(S))),'fro');
        if a>delta
            k = k;
            break;
        end
    end
    S_trun = S(1:k,1:k);
    U_trun = U(:,1:k);
    V_trun = V(:,1:k);
    C = U_trun*S_trun*V_trun';
    r(j) = rank(C);
    TTt.cores{j,1} = reshape(V_trun',[r(j),n(j),r(j+1)]);
    TTt.dim(j,:) = [r(j),n(j),r(j+1)];
    C = U_trun*S_trun;
end
a = TTn.dim(1,:);
TTt.cores{1,1} = permute(reshape(TTn.cores{1,1},[a(1)*a(2),a(3)])*C,[3,1,2]);
TTt.dim(1,:) = [1,a(1)*a(2),size(C,2)];
%% checked on 04/05/2021
end


% TTori = TT;
% [d,~] = size(TT);
% norm1 = innerprodTT(TT,TT);
% delta = eps/sqrt(d-1)*sqrt(norm1);
% 
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
% 
% % the norm of the tensor train is now in the last tensor train matrix.
% 
% % right to left SVD-decomposition procedure
% r = zeros(1,d+1);r(1) = 1;r(end) = 1;
% n = zeros(d,1);
% a = size(TTn{d,1});
% n(d) = a(2);
% [U,S,V] = svd(TTn{d,1},'econ');
% 
% % delta-truncate SVD, take frobinius norm of singular values, E is what
% % you throw away.
% for k = min(size(S)):-1:1
%     a = norm(S(k:min(size(S)),k:min(size(S))),'fro');
%     if a>delta
%         k = k;
%         break;
%     end
% end
% S_trun = S(1:k,1:k);
% U_trun = U(:,1:k);
% V_trun = V(:,1:k);
% 
% C = U_trun*S_trun*V_trun';
% r(d) = rank(C);
% TTt{d,1} = reshape(V_trun',[r(d),n(d),r(d+1)]); 
% C = U_trun*S_trun;
% for j = d-1:-1:2
%     a = size(TTn{j,1});
%     n(j) = a(2);
%     b = size(C);
%     TTe = reshape(reshape(TTn{j,1},[a(1)*a(2),a(3)])*C,[a(1),a(2)*b(2)]);
%     [U,S,V] = svd(TTe,'econ');
% 
%     % delta-truncate SVD, take frobinius norm of singular values, E is what
%     % you throw away.
%     for k = min(size(S)):-1:1
%         a = norm(S(k:min(size(S)),k:min(size(S))),'fro');
%         if a>delta
%             k = k;
%             break;
%         end
%     end
%     S_trun = S(1:k,1:k);
%     U_trun = U(:,1:k);
%     V_trun = V(:,1:k);
%     C = U_trun*S_trun*V_trun';
%     r(j) = rank(C);
%     TTt{j,1} = reshape(V_trun',[r(j),n(j),r(j+1)]); 
%     C = U_trun*S_trun;
% end
% a = size(TTn{1,1});
% TTt{1,1} = permute(reshape(TTn{1,1},[a(1)*a(2),a(3)])*C,[3,1,2]);