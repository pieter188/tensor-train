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


%%
% a = size(P);
% v2 = zeros(size(P)); % orthogonal basis
% q2 = zeros(size(P)); % orthonormal basis
% for i = 1:a(2)
%     v2(:,i) = P(:,i);
% end
% for i = 1:a(2)
%     q2(:,i) = v2(:,i);
%     for j = i+1:a(2)
%         v2(:,j) = v2(:,j) - ((q2(:,i)'*v2(:,j))/(q2(:,i)'*q2(:,i)))*q2(:,i);
%     end
%     q2(:,i) = v2(:,i)/sqrt(v2(:,i)'*v2(:,i));
% end
% r2 = q2'*P;


%% Modified Gram-Schmidt algorithm in Tensor Train format for 20 by 20 matrix
% Row = [3,3,4,2,2];
% Col = [3,3,4,2,2];
% Pten = Mat2Tensor(P,Row,Col);
% eps = 1e-6;
% PTT = Tensor2TT_SVD_eps(Pten,eps);
% PTTm = TT2TTm(PTT,Row,Col);

for alpha = 1:10
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
        rank_threshold = 90;
        if max(V_tp{1,j}.dim(:,1)) > rank_threshold
            eps = 1*10^(1-alpha);
            V_test = V_tp{1,j};
            [V_tp{1,j},~] = roundTT2(V_tp{1,j},eps);
        end
    end   
end
% %% obtain data
% % initialize table
% T(1,:) = {'Qi','TT-rank Qi','Qj','TT-rank Qj','innerproduct; Qi*Qj',...
%     'epsilon','TT-cores','Row Partitioning','Column Partitioning',...
%     'Rank Threshold','innerproduct matrix case: qi*qj'};
% h = 2;
% for i = 1:prod(Col)
%     for j = i:prod(Col)
%         % obtain ranks
%         r_i(1:d) = Q2{1,i}.dim(:,1); r_i(d+1) = 1;
%         r_j(1:d) = Q2{1,j}.dim(:,1); r_j(d+1) = 1;
%         % obtain innerproduct
%         inner_ij = innerprodTT(Q2{1,i},Q2{1,j});
%         inner_mat = q1(:,i)'*q1(:,j);
%         T(h,:) = {i, r_i , j, r_j, inner_ij, eps, d, Row,...
%             Col, rank_threshold,inner_mat}; 
%         h = h + 1;
%     end
% end
% 
% %% import data to excel spreadsheet
% filename = 'Testingdata_144x144_150_Real_ortho_test.xlsx'; %input('Please define the name of the spreadsheet. This can also be an existing spreadsheet. (Must be a .xlsx file)');
% sheetnumber = alpha;
% writecell(T,filename,'Sheet',sheetnumber,'range','A1')
end