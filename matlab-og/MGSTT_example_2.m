%% Modified Gram-Schmidt algorithm for 1000 by 1000 matrix
% Pieter van Klaveren
% 01-04-2021
clear all; clc; close all

A = rand(2^9,2^9)*5;
a = size(A);
v1 = zeros(size(A)); % orthogonal basis
q1 = zeros(size(A)); % orthonormal basis
for i = 1:a(2)
    v1(:,i) = A(:,i);
end
for i = 1:a(2)
    q1(:,i) = v1(:,i)/sqrt(v1(:,i)'*v1(:,i));
    for j = i+1:a(2)
        v1(:,j) = v1(:,j) - (q1(:,i)'*v1(:,j))*q1(:,i);
    end
end

%% Modified Gram-Schmidt algorithm in Tensor Train format for 1000 by 1000 matrix
Row = 2^3*ones(1,3);
Col = 2^3*ones(1,3);
Aten = Mat2Tensor(A,Row,Col);
eps = 1e-16;
ATT = Tensor2TT_SVD_eps(Aten,eps);
ATTm = TT2TTm(ATT,Row,Col);
clear Aten ATT
for alpha = 1:10
    for i = 1:prod(Col)
        % take ith vector out of TTm of A
        % function to determine the code of indexing.
        sizeATTm = size(ATTm);
        e{1,1} = zeros(Col(1),1);
        e{2,1} = zeros(Col(2),1);
        e{3,1} = zeros(Col(3),1);
        % determine e3
        for s = 1:Col(3)
            if s-1 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= s
                e3 = s;
                break;
            end
        end
        e{3,1}(e3,1) = 1;
        % determine e2
        if (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 1
            e2 = 1;
        else 
            e2 = 2;
        end
        e{2,1}(e2,1) = 1;

        % determine e1
        if i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 1
            e1 = 1;
        else 
            e1 = 2;
        end
        e{1,1}(e1,1) = 1;

        % column pick from ATTm via e1,e2,e3
        [d,~] = size(ATTm);
        for k = 1:d
            a = size(ATTm{k,1});
            if k == d
               a = [a,1]; 
            end
            V{1,i}{k,1} = reshape(permute(ATTm{k,1},[1,2,4,3]),[numel(ATTm{k,1})/a(3),a(3)])*e{k,1};
            V{1,i}{k,1} = reshape(V{1,i}{k,1},[a(1),a(2),a(4)]);
        end
    end   
        % this for-loop can be divided into three segments: norm
        % computation of jth column and of jth with ith column, multiplying
        % jth vector with this norm ratio and substraction.
        % norm of jth column, and ith*jth column
    for i = 1:prod(Col)
        Q1{1,i} = V{1,i};
        Q1{1,i} = sitek(Q1{1,i},1);
        normQ1 = norm(reshape(Q1{1,i}{1,1},[numel(Q1{1,i}{1,1}),1]),'fro');
        Q1{1,i}{1,1} = Q1{1,i}{1,1}/normQ1;
        for j = i+1:prod(Col)
            [normij,~] = innerprodTT(Q1{1,i},V{1,j});
            % subtract jth column from ith column, this can be done via
            % addition of v(:,i) and -v(:,j)
            QQ = Q1{1,i};
            QQ{1,1} = normij*QQ{1,1};
            V{1,j} = subTT(V{1,j},QQ);
            % make function that outputs the ranks
            rc(1) = 1;rc(d+1)=1;
            for p = 1:d-1
                c = size(V{1,j}{p,1});
                rc(p+1) = c(3);
            end
            rank_threshold = 130;
            if max(rc) > rank_threshold
%                 r = [1,8,8,1];
%                 [V{1,j},rt] = roundTTr(V{1,j},r);
                eps = 1*10^(-alpha+1);
                [V{1,j},r] = roundTT(V{1,j},eps);
            end
        end    
    end
    %% obtain data
    % initialize table
    T(1,:) = {'Qi','TT-rank Qi','Qj','TT-rank Qj','innerproduct; Qi*Qj','epsilon','TT-cores','Row Partitioning','Column Partitioning','Rank Threshold'};
    h = 2;
    for i = 1:prod(Col)
        for j = i:prod(Col)
            % obtain ranks
            r_i(1) = 1; r_i(d+1)=1;
            r_j(1) = 1; r_j(d+1)=1;
            for p = 1:d-1
                c = size(Q1{1,i}{p,1});
                if length(c) == 2
                    c = [c,1];
                end
                r_i(p+1) = c(3);
                c = size(Q1{1,j}{p,1});
                if length(c) == 2
                    c = [c,1];
                end
                r_j(p+1) = c(3);
            end
            % obtain innerproduct
            [inner_ij,~] = innerprodTT(Q1{1,i},Q1{1,j});
            T(h,:) = {i, r_i , j, r_j, inner_ij, eps, d, Row, Col, rank_threshold}; 
            h = h + 1;
        end
    end

    %% import data to excel spreadsheet
    filename = 'Testingdata_100x100.xlsx'; %input('Please define the name of the spreadsheet. This can also be an existing spreadsheet. (Must be a .xlsx file)');
    sheetnumber = alpha;
    writecell(T,filename,'Sheet',sheetnumber,'range','A1')
end