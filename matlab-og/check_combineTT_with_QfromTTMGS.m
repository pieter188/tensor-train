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


%% Modified Gram-Schmidt algorithm in Tensor Train format for 20 by 20 matrix
% Row = [3,3,4,2,2];
% Col = [3,3,4,2,2];
% Pten = Mat2Tensor(P,Row,Col);
% eps = 1e-6;
% PTT = Tensor2TT_SVD_eps(Pten,eps);
% PTTm = TT2TTm(PTT,Row,Col);

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
%%
for alpha = 5:5
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
        eps(alpha) = 1*10^(-alpha+1);
%         Vtest = Vf;
        [Vf,~] = roundTT2(Vf,eps(alpha)); % kijk hier naar
%         aaa = aaa+1;
%         Vtest1 = TT2Tensor(Vtest);
%         Vtest1 = Tensor2Mat(Vtest1,[3 3 4 2 2],[3 3 4 2 2]);
%         Vf1 = TT2Tensor(Vf);
%         Vf1 = Tensor2Mat(Vf1,Row,Col);
%         norm((Vtest1-Vf1),'fro')/norm(Vtest1,'fro')
    end
end
Vf = TT2TTm(Vf,Row1,Col1);
% Vf = combineTT(V_tp,Col);
% Vf = TTm2TT(Vf);
% Vf = TT2Tensor(Vf);
% Vf = Tensor2Mat(Vf,Row,Col);
%---------------------------------------------------------------------
% errorQ350_norm(alpha) = norm(eye(144)-(Vf'*Vf),'fro')/norm(eye(144),'fro');

end

%% figure
% close all
% load('errorQ5_norm')
% load('errorQ20_norm')
% load('errorQ50_norm')
% load('errorQ90_norm')
% load('errorQ150_norm')
% load('errorQ200_norm')
% load('errorQ250_norm')
% % load('errorQ300_norm')
% figure
% loglog(eps(1:7),errorQ50_norm(1:7),'k-*','LineWidth',1.5)
% hold on
% loglog(eps(1:7),errorQ150_norm(1:7),'b-d','LineWidth',1.5)
% loglog(eps(1:7),errorQ250_norm(1:7),'r-o','LineWidth',1.5)
% % loglog(eps,errorQ300_norm,'k-->','LineWidth',1.5)
% grid on
% xlabel('Epsilon')
% ylabel('Orthogonality error')
% title('Orthogonality error of combineTT function')
% legend('Rank threshold: 50',...
%     'Rank threshold: 150',...
%     'Rank threshold: 250','location','northwest');