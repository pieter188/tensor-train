%% Testing the validity of the orthogonal error
% Pieter van Klaveren
% 01-04-2021
clear all; clc; close all
res = [9,16];
sigma_p = 1;
P = eye(prod(res))*sigma_p^2;
load('State_covmatrix_tp_test.mat')
PTTm = LL{1,53};
[P,Row,Col] = TTm2TT(PTTm);
P = TT2Tensor(P);
P = Tensor2Mat(P,Row,Col);
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
Q = q1(:,1);
for i = 2:a(2)
    Q = [Q,q1(:,i)];
end
R = Q'*P;

%% Modified Gram-Schmidt algorithm in Tensor Train format

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
        rank_threshold = 30;
        if max(V_tp{1,j}.dim(1,:)) > rank_threshold
            eps = 1*10^(-16);
            [V_tp{1,j},~] = roundTT2(V_tp{1,j},eps);
        end
    end   
end
for alpha = 1:10
clear Vm
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
for i = 2:prod(Cdim)
    Vf = addTT(Vf,Vm{1,i});
%     Vf = sitek(Vf,1);
    rank_threshold = 30;
    if max(Vf.dim(:,1)) > rank_threshold
        eps = 1*10^(-alpha+1);
        [Vf,~] = roundTT2(Vf,eps); % kijk hier naar
    end
end
Vf = TT2TTm(Vf,Row1,Col1);
Q = Vf;
%---------------------------------------------------------------------


%     Q = combineTT(Q2,Col);
RTTm = MatMatTT(TransposeTTm(Q),PTTm);
[RTTm,RowR,ColR] = TTm2TT(RTTm);
RTTm = TT2Tensor(RTTm);
RTTm = Tensor2Mat(RTTm,RowR,ColR);
errorR = norm(R-RTTm,'fro')/norm(R,'fro');
%% obtain data
% initialize table
T(1,:) = {'Qi','TT-rank Qi','Qj','TT-rank Qj','innerproduct; Qi*Qj',...
    'epsilon','TT-cores','Row Partitioning','Column Partitioning',...
    'Rank Threshold','innerproduct matrix case: qi*qj','error in R'};
h = 2;
for i = 1:prod(Col)
    for j = i:prod(Col)
        % obtain ranks
        r_i(1:d) = Q2{1,i}.dim(:,1); r_i(d+1) = 1;
        r_j(1:d) = Q2{1,j}.dim(:,1); r_j(d+1) = 1;
        % obtain innerproduct
        inner_ij = innerprodTT(Q2{1,i},Q2{1,j});
        inner_mat = q1(:,i)'*q1(:,j);
        T(h,:) = {i, r_i , j, r_j, inner_ij, eps, d, Row,...
            Col, rank_threshold,inner_mat,errorR}; 
        h = h + 1;
    end
end


%% import data to excel spreadsheet
filename = 'Testingdata_144x144_30_Real_recR_combineTT_allinfokept.xlsx'; %input('Please define the name of the spreadsheet. This can also be an existing spreadsheet. (Must be a .xlsx file)');
sheetnumber = alpha;
writecell(T,filename,'Sheet',sheetnumber,'range','A1')
end
