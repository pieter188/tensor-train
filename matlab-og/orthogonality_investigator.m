%% Orthogonality checker
% Pieter van Klaveren
% 14/05/2021
clear all; close all; clc
% This function checks whether a TT-representation of a certain orthogonal
% matrix is still orthogonal, for which ranks or truncation parameter is the
% matrix still the orthogonal.

%% create orthogonal matrix
% load('State_covmatrix_testpurpose.mat')
% load('State_covmatrix_tp_test.mat')
% L = LLL{1,3};
% [L,Row,Col] = TTm2TT(L);
% L = TT2Tensor(L);
% L = Tensor2Mat(L,Row,Col);
% if rank(L) == size(L,1)
%     display("L is full rank")
% else
%     display("L is not full rank")
% end

L = rand(20);
Row = [2 2 5];
Col = [2 2 5];

%% create TTm representation of the matrix and check its orthogonality for eps
for alpha = 1:20
    % Q is the orthogonal basis for L
    Q = orth(L);
    Q = Mat2Tensor(Q,Row,Col);
    eps(alpha) = 1*10^(1-alpha);
    Q = Tensor2TT_SVD_eps(Q,eps(alpha));
    Q = TT2TTm(Q,Row,Col);

    I = MatMatTT(TransposeTTm(Q),Q);    % ranks are very high
    I = TTm2TT(I);
    I = TT2Tensor(I);
    I = Tensor2Mat(I,Row,Col);
    I = I-eye(prod(Row));

    I = tril(I);
    I = reshape(I,[numel(I),1]);
    I(I==0) = [];
    Imean(alpha) = mean(abs(I));
    Imax(alpha) = max(abs(I));
end

%% make a boxplot of the errors encountered while checking
figure
loglog(eps,Imean,'k-')
hold on
loglog(eps,Imax,'r-')
grid on
ylabel('error')
xlabel('eps')
legend('mean error','max error')

%% create TTm representation of the matrix and check its orthogonality for rank
for alpha = 1:200
    % Q is the orthogonal basis for L
    Q = orth(L);
    Q = Mat2Tensor(Q,Row,Col);
    eps = 1e-20;
    Q = Tensor2TT_SVD_eps(Q,eps);
    Q = roundTTr(Q,[1,min(Q.dim(2,1),alpha),min(Q.dim(3,1),alpha),...
        1]);
    r_max(alpha) = max([1,min(Q.dim(2,1),alpha),min(Q.dim(3,1),alpha),...
        1]);
    Q = TT2TTm(Q,Row,Col);
    
    I = MatMatTT(TransposeTTm(Q),Q);    % ranks are very high
    I = TTm2TT(I);
    I = TT2Tensor(I);
    I = Tensor2Mat(I,Row,Col);
    error1(alpha) = norm(I-eye(size(L,1)),'fro')/norm(eye(size(L,1)),'fro');
    
    
%     I = I-eye(prod(Row));
% 
%     I = tril(I);
%     I = reshape(I,[numel(I),1]);
%     I(I==0) = [];
%     Imeanr(alpha) = mean(abs(I));
%     Imaxr(alpha) = max(abs(I));
    
    
end

%% make a boxplot of the errors encountered while checking
figure
semilogy(r_max,error1,'k-')
% hold on
% semilogy(r_max,Imaxr,'r-')
grid on
ylabel('error')
xlabel('maximum rank')
legend('mean error','max error')

%% random matrix
L = rand(size(L,1))+rand(size(L,1))+rand(size(L,1))/0.5;
for alpha = 1:200
    Q = L;
    Q = Mat2Tensor(Q,Row,Col);
    eps = 1e-20;
    Q = Tensor2TT_SVD_eps(Q,eps);
    Q = roundTTr(Q,[1,min(Q.dim(2,1),alpha),min(Q.dim(3,1),alpha),...
        1]);
    r_max1(alpha) = max([1,min(Q.dim(2,1),alpha),min(Q.dim(3,1),alpha),...
        1]);
    
    Q = TT2Tensor(Q);
    Q = Tensor2Mat(Q,Row,Col);
    
    error(alpha) = norm(L-Q,'fro')/norm(L,'fro');
end

figure
semilogy(r_max1,error,'k-')
grid on
ylabel('error')
xlabel('maximum rank')
legend('error')
title('error normal matrix with rank truncation')





%% rank
clear all; close all
load('orthogonalmat_testpurpose.mat')
load('State_covmatrix_testpurpose.mat')
load('State_covmatrix_tp_test.mat')
L = LL{1,53};
Lmat = TTm2TT(L);
Lmat = TT2Tensor(Lmat);
Lmat = Tensor2Mat(Lmat,[3 3 4 2 4],[3 3 4 2 2]);
Qmat = orth(Lmat);
Q = Mat2Tensor(Qmat,[3 3 4 2 4],[3 3 4 2 2]);
eps = 1*10^(-20);
Q = Tensor2TT_SVD_eps(Q,eps);
Q = TT2TTm(Q,[3 3 4 2 4],[3 3 4 2 2]);

% Q = QQQ{1,5};
% Qmat = TTm2TT(Q);
% Qmat = TT2Tensor(Qmat);
% Qmat = Tensor2Mat(Qmat,[3 3 4 2 4],[3 3 4 2 2]);

% Qmat = rand(288,144);
for alpha = 1:100
    Qtrun = TTm2TT(Q);
    Qtrun = roundTTr(Qtrun,[1,min(Qtrun.dim(2,1),alpha),min(Qtrun.dim(3,1),alpha),...
            min(Qtrun.dim(4,1),alpha), min(Qtrun.dim(5,1),alpha), 1]);
    r_max1(alpha) = max([1,min(Q.dim(2,1),alpha),min(Q.dim(3,1),alpha),...
            min(Q.dim(4,1),alpha), min(Q.dim(5,1),alpha), 1]);
    Qtrun = TT2Tensor(Qtrun);
    Qtrun = Tensor2Mat(Qtrun,[3 3 4 2 4],[3 3 4 2 2]);
    error(alpha) = norm(Qmat-Qtrun,'fro')/norm(Qmat,'fro');
    Ortho_error(alpha) = norm(eye(144)-Qtrun'*Qtrun,'fro')/norm(eye(144),'fro');
end
Ortho_error_notrunc = norm(eye(144)-Qmat'*Qmat,'fro')/norm(eye(144),'fro');

figure
semilogy(r_max1,error,'k-d')
grid on
title('Relative error on truncated orthogonal matrix Q relative to original Q')
xlabel('maximum rank')
ylabel('relative error')
legend('relative error of truncated Q to original Q')

figure
semilogy(r_max1,Ortho_error,'k-d')
hold on
semilogy(max(r_max1),Ortho_error_notrunc,'r*')
grid on
title('Relative error on truncated orthogonal matrix Q^T*Q relative to I')
xlabel('maximum rank')
ylabel('relative error')
legend('relative error to I with truncation','relative error to I no truncation')

%% epsilon
clear all; close all
load('orthogonalmat_testpurpose.mat')
% load('State_covmatrix_testpurpose.mat')
load('State_covmatrix_tp_test.mat')
L = LL{1,53};
Lmat = TTm2TT(L);
Lmat = TT2Tensor(Lmat);
Lmat = Tensor2Mat(Lmat,[3 3 4 2 4],[3 3 4 2 2]);
Qmat = orth(Lmat);
Q = Mat2Tensor(Qmat,[3 3 4 2 4],[3 3 4 2 2]);
eps = 1*10^(-20);
Q = Tensor2TT_SVD_eps(Q,eps);
Q = TT2TTm(Q,[3 3 4 2 4],[3 3 4 2 2]);

% Q = QQQ{1,5};
% Qmat = TTm2TT(Q);
% Qmat = TT2Tensor(Qmat);
% Qmat = Tensor2Mat(Qmat,[3 3 4 2 4],[3 3 4 2 2]);

for alpha = 1:20
    Qtrun = TTm2TT(Q);
    eps(alpha) = 1*10^(1-alpha);
    Qtrun = roundTT2(Qtrun,eps(alpha));
    Qtrun = TT2Tensor(Qtrun);
    Qtrun = Tensor2Mat(Qtrun,[3 3 4 2 4],[3 3 4 2 2]);
    error(alpha) = norm(Qmat-Qtrun,'fro')/norm(Qmat,'fro');
    Ortho_error(alpha) = norm(eye(144)-Qtrun'*Qtrun,'fro')/norm(eye(144),'fro');
end
Ortho_error_notrunc = norm(eye(144)-Qmat'*Qmat,'fro')/norm(eye(144),'fro');

figure
loglog(eps,error,'k-d')
grid on
title('Relative error on truncated orthogonal matrix Q relative to original Q')
xlabel('epsilon')
ylabel('relative error')
legend('relative error of truncated Q to original Q')

figure
loglog(eps,Ortho_error,'k-d')
hold on
loglog(10^(-19),Ortho_error_notrunc,'r*')
grid on
title('Relative error on truncated orthogonal matrix Q^T*Q relative to I')
xlabel('epsilon')
ylabel('relative error')
legend('relative error to I with truncation','relative error to I no truncation')

%% Modified Gram-Schmidt algorithm for 20 by 20 matrix
% Pieter van Klaveren
% 01-04-2021
clear all; clc; close all
res = [9,16];
% sigma_p = 1;
% P = eye(prod(res))*sigma_p^2;
load('State_covmatrix_testpurpose.mat')
PTTm = LLL{1,53};
[P,Row,Col] = TTm2TT(PTTm);
P = TT2Tensor(P);
P = Tensor2Mat(P,Row,Col);

% P = (rand(144)+rand(144)+randn(144));
a = size(P);
v1 = zeros(size(P)); % orthogonal basis
q1 = zeros(size(P)); % orthonormal basis
for i = 1:a(2)
    v1(:,i) = P(:,i);
end
for i = 1:a(2)
    q1(:,i) = v1(:,i)/sqrt(v1(:,i)'*v1(:,i));
    for j = i+1:a(2)
        v1(:,j) = v1(:,j) - (q1(:,i)'*v1(:,j))*q1(:,i);
    end
end
r1 = q1'*P;

for i = 1:prod(Col)
    % take ith vector out of TTm of A
    % function to determine the code of indexing.
    e = multi_index(Col,i);
    % column pick from ATTm via e1,e2,e3
    d = size(PTTm.dim,1);
    for m = 1:d
        V.cores{m,1} = reshape(permute(PTTm.cores{m,1},[1,2,4,3]),...
            [numel(PTTm.cores{m,1})/PTTm.dim(m,3),PTTm.dim(m,3)])*e{m,1};
        V.cores{m,1} = reshape(V.cores{m,1},...
            [PTTm.dim(m,1),PTTm.dim(m,2),PTTm.dim(m,4)]);
        V.dim(m,:) = [PTTm.dim(m,1),PTTm.dim(m,2),PTTm.dim(m,4)];
    end
    V_tp{1,i} = V;
end     
% make vectors orthogonal to each other
for i = 1:prod(Col)
    % normalize Q via site-1 orthogonalization and division by first
    % core norm
    Q2{1,i} = V_tp{1,i};
    Q2{1,i} = sitek(Q2{1,i},1);
    normQ1 = norm(reshape(Q2{1,i}.cores{1,1},[numel(Q2{1,i}.cores...
        {1,1}),1]),'fro');
    Q2{1,i}.cores{1,1} = Q2{1,i}.cores{1,1}/normQ1;
    for j = i+1:prod(Col)
        normij = innerprodTT(Q2{1,i},V_tp{1,j});
        QQ = Q2{1,i};
        QQ.cores{1,1} = normij*QQ.cores{1,1};
        V_tp{1,j} = subTT(V_tp{1,j},QQ);
%             V_tp{1,j} = sitek(V_tp{1,j},1);
        % check rank and truncate
        rank_threshold = 5;
        if max(V_tp{1,j}.dim(:,1)) > rank_threshold
            eps = 1*10^(-6);
            V_test = V_tp{1,j};
            [V_tp{1,j},~] = roundTT2(V_tp{1,j},eps);
        end
    end   
end
% output here: set of orthonormal vectors stored in struct Q2.
% --------------------------------------------------------------------
clear Vm V
V = Q2;
Cdim = [3 3 4 2 2];
%---------------------------------------------------------------------
d = size(V{1,1}.dim,1);
for i = 1:prod(Cdim)
    E = multi_index(Cdim,i);
    for k = 1:d
        a = V{1,i}.dim(k,:);
        b = size(E{k,1});
        Vm{1,i}.cores{k,1} = reshape(V{1,i}.cores{k,1},[a(1),a(2),1,a(3)]);
        Vm{1,i}.cores{k,1} = reshape(permute(Vm{1,i}.cores{k,1},[1,2,4,3]),...
            [a(1)*a(2)*a(3),1])*E{k,1}';
        Vm{1,i}.cores{k,1} = reshape(Vm{1,i}.cores{k,1},[a(1),a(2),a(3),b(1)]);
        Vm{1,i}.cores{k,1} = permute(Vm{1,i}.cores{k,1},[1,2,4,3]);
        Vm{1,i}.dim(k,:) = [a(1),a(2),b(1),a(3)];
    end
    [Vm{1,i},Row1,Col1] = TTm2TT(Vm{1,i});
end
% sum the vectors up in TT format
Vf = Vm{1,1};
aaa = 0;
for i = 2:prod(Cdim)
    Vf = addTT(Vf,Vm{1,i});
    rank_threshold = 150;
    if max(Vf.dim(:,1)) > rank_threshold
        eps = 1*10^(-4);
%         Vtest = Vf;
        [Vf,~] = roundTT2(Vf,eps); % kijk hier naar
%         aaa = aaa+1;
%         Vtest1 = TT2Tensor(Vtest);
%         Vtest1 = Tensor2Mat(Vtest1,[3 3 4 2 2],[3 3 4 2 2]);
%         Vf1 = TT2Tensor(Vf);
%         Vf1 = Tensor2Mat(Vf1,Row,Col);
%         norm((Vtest1-Vf1),'fro')/norm(Vtest1,'fro')
    end
end
Vf = TT2TTm(Vf,Row1,Col1);
Q = Vf;
Qmat = TTm2TT(Q);
Qmat = TT2Tensor(Qmat);
Qmat = Tensor2Mat(Qmat,[3,3,4,2,2],[3,3,4,2,2]);
% rank
close all
clear Qtrun eps error Ortho_error Ortho_error_notrunc 
% for alpha = 1:100
%     Qtrun = TTm2TT(Q);
%     Qtrun = roundTTr(Qtrun,[1,min(Qtrun.dim(2,1),alpha),min(Qtrun.dim(3,1),alpha),...
%             min(Qtrun.dim(4,1),alpha), min(Qtrun.dim(5,1),alpha), 1]);
%     r_max1(alpha) = max([1,min(Q.dim(2,1),alpha),min(Q.dim(3,1),alpha),...
%             min(Q.dim(4,1),alpha), min(Q.dim(5,1),alpha), 1]);
%     Qtrun = TT2Tensor(Qtrun);
%     Qtrun = Tensor2Mat(Qtrun,[3 3 4 2 2],[3 3 4 2 2]);
%     error(alpha) = norm(Qmat-Qtrun,'fro')/norm(Qmat,'fro');
%     Ortho_error(alpha) = norm(eye(144)-Qtrun'*Qtrun,'fro')/norm(eye(144),'fro');
% end
Qtrun = TTm2TT(Q);
Qtrun = roundTTr(Qtrun,[1,9,80,16,4,1]);
Qtrun1 = Qtrun;
r_max1 = max(Qtrun.dim(:,1));
Qtrun = TT2Tensor(Qtrun);
Qtrun = Tensor2Mat(Qtrun,[3 3 4 2 2],[3 3 4 2 2]);
error = norm(Qmat-Qtrun,'fro')/norm(Qmat,'fro');

% figure
% semilogy(r_max1,error,'k-d')
% grid on
% title('Relative error on truncated orthogonal matrix Q relative to original Q')
% xlabel('maximum rank')
% ylabel('relative error')
% legend('relative error of truncated Q to original Q')

% figure
% semilogy(r_max1,Ortho_error,'k-d')
% hold on
% semilogy(max(r_max1),Ortho_error_notrunc,'r*')
% grid on
% title('Relative error on truncated orthogonal matrix Q^T*Q relative to I')
% xlabel('maximum rank')
% ylabel('relative error')
% legend('relative error to I with truncation','relative error to I no truncation')

%% epsilon
close all
clear Qtrun eps error Ortho_error Ortho_error_notrunc 
for alpha = 1:20
    Qtrun = TTm2TT(Q);
    eps(alpha) = 1*10^(1-alpha);
    Qtrun = roundTT2(Qtrun,eps(alpha));
    Qtrun = TT2Tensor(Qtrun);
    Qtrun = Tensor2Mat(Qtrun,[3 3 4 2 2],[3 3 4 2 2]);
    error(alpha) = norm(Qmat-Qtrun,'fro')/norm(Qmat,'fro');
    Ortho_error(alpha) = norm(eye(144)-Qtrun'*Qtrun,'fro')/norm(eye(144),'fro');
end
Ortho_error_notrunc = norm(eye(144)-Qmat'*Qmat,'fro')/norm(eye(144),'fro');

figure
loglog(eps,error,'k-d','LineWidth',1.5)
grid on
title('Relative error on truncated orthogonal matrix Q relative to original Q')
xlabel('epsilon')
ylabel('relative error')
legend('relative error of truncated Q to original Q','location','northwest')

figure
loglog(eps,Ortho_error,'k-d','LineWidth',1.5)
hold on
loglog(10^(-19),Ortho_error_notrunc,'r*','LineWidth',1.5)
grid on
title('Relative error on truncated orthogonal matrix Q^T*Q relative to I')
xlabel('epsilon')
ylabel('relative error')
legend('relative error to I with truncation','relative error to I no truncation','location','northwest')


















