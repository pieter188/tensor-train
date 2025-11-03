%% Testing the validity of the orthogonal error
% Pieter van Klaveren
% 01-04-2021
clear all; clc; close all
res = [9,16];
sigma_p = 1;
P = eye(prod(res))*sigma_p^2;
load('State_covmatrix_tp_test.mat')
load('State_covmatrix_testpurpose.mat')
PTTm = LLL{1,53};
[P,Row,Col] = TTm2TT(PTTm);
P = TT2Tensor(P);
P = Tensor2Mat(P,Row,Col);
Q = orth(P);
QQ = Q;
Q = Mat2Tensor(Q,Row,Col);
eps = 1e-20;
Q = Tensor2TT_SVD_eps(Q,eps);
Q = TT2TTm(Q,Row,Col);
%%
for i = 1:prod(Col)
    % take ith vector out of TTm of A
    % function to determine the code of indexing.
    e = multi_index(Col,i);
    % column pick from ATTm via e1,e2,e3
    d = size(Q.dim,1);
    for m = 1:d
        V.cores{m,1} = reshape(permute(Q.cores{m,1},[1,2,4,3]),...
            [numel(Q.cores{m,1})/Q.dim(m,3),Q.dim(m,3)])*e{m,1};
        V.cores{m,1} = reshape(V.cores{m,1},...
            [Q.dim(m,1),Q.dim(m,2),Q.dim(m,4)]);
        V.dim(m,:) = [Q.dim(m,1),Q.dim(m,2),Q.dim(m,4)];
    end
    V_tp{1,i} = V;
end     

for alpha = 1:20
clear Vm V
V = V_tp;
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
for i = 2:prod(Cdim)
    Vf = addTT(Vf,Vm{1,i});
%     Vf = sitek(Vf,1);
    rank_threshold = 90;
    if max(Vf.dim(:,1)) > rank_threshold
        eps(alpha) = 1*10^(-alpha+1);
        Vtest = Vf;
        [Vf,~] = roundTT2(Vf,eps(alpha)); % kijk hier naar
    end
end
Vf = TT2TTm(Vf,Row1,Col1);
% Vf = combineTT(V_tp,Col);
Vf = TTm2TT(Vf);
Vf = TT2Tensor(Vf);
Vf = Tensor2Mat(Vf,Row,Col);
%---------------------------------------------------------------------
errorQ(alpha) = norm(QQ-Vf,'fro')/norm(QQ,'fro');

end

%% figure
close all
load('errorQ5')
load('errorQ10')
load('errorQ20')
load('errorQ30')
load('errorQ50')
load('errorQ70')
load('errorQ90')
load('errorQ120')
figure
loglog(eps,errorQ5,'k-*','LineWidth',1.5)
hold on
loglog(eps,errorQ10,'b-^','LineWidth',1.5)
loglog(eps,errorQ20,'g-+','LineWidth',1.5)
loglog(eps,errorQ30,'y-d','LineWidth',1.5)
loglog(eps,errorQ50,'m-h','LineWidth',1.5)
loglog(eps,errorQ70,'c-<','LineWidth',1.5)
loglog(eps,errorQ90,'r-o','LineWidth',1.5)
loglog(eps,errorQ120,'k-->','LineWidth',1.5)
grid on
xlabel('epsilon')
ylabel('relative error of orthogonal matrix')
title('relative error of orthogonal matrix for combineTT function for varying rank threshold and epsilon')
legend('Rank threshold: 5','Rank threshold: 10','Rank threshold: 20',...
    'Rank threshold: 30','Rank threshold: 50','Rank threshold: 70',...
    'Rank threshold: 90','Rank threshold: 120');