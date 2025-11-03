%% Gram-Schmidt algorithm for 20 by 20 matrix
% Pieter van Klaveren
% 01-04-2021
clear all; clc; close all
for xxx = 1:1
% A = rand(20,20)*5;
% a = size(A);
% v = zeros(size(A)); % orthogonal basis
% q = zeros(size(A)); % orthonormal basis
% for i = 1:a(2)
%     v(:,i) = A(:,i);
%     for j = 1:i-1
%         v(:,i) = v(:,i) - v(:,j)'*v(:,i)/(v(:,j)'*v(:,j))*v(:,j);
%     end    
%     q(:,i) = v(:,i)/norm(v(:,i));
% end
% 
% %% Gram-Schmidt algorithm in Tensor Train format for 20 by 20 matrix
% Row = [2,2,5];
% Col = [2,2,5];
% Aten = Mat2Tensor(A,Row,Col);
% eps = 0.01;
% ATT = Tensor2TT_SVD_eps(Aten,eps);
% ATTm = TT2TTm(ATT,Row,Col);
% clear Aten ATT
% 
% for i = 1:prod(Col)
%     % take ith vector out of TTm of A
%     % function to determine the code of indexing.
%     sizeATTm = size(ATTm);
%     e{1,1} = zeros(Col(1),1);
%     e{2,1} = zeros(Col(2),1);
%     e{3,1} = zeros(Col(3),1);
%     % determine e3
%     if i/Col(1)/Col(2) <= 1
%         e3 = 1;
%     elseif 1 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 2
%         e3 = 2;
%     elseif 2 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 3
%         e3 = 3;
%     elseif 3 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 4
%         e3 = 4;
%     else 
%         e3 = 5;
%     end
%     e{3,1}(e3,1) = 1;
%     
%     % determine e2
%     if (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 1
%         e2 = 1;
%     else 
%         e2 = 2;
%     end
%     e{2,1}(e2,1) = 1;
%     
%     % determine e1
%     if i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 1
%         e1 = 1;
%     else 
%         e1 = 2;
%     end
%     e{1,1}(e1,1) = 1;
%     
%     % column pick from ATTm via e1,e2,e3
%     [d,~] = size(ATTm);
%     for k = 1:d
%         a = size(ATTm{k,1});
%         if k == d
%            a = [a,1]; 
%         end
%         V{1,i}{k,1} = reshape(permute(ATTm{k,1},[1,2,4,3]),[numel(ATTm{k,1})/a(3),a(3)])*e{k,1};
%         V{1,i}{k,1} = reshape(V{1,i}{k,1},[a(1),a(2),a(4)]);
%     end
%     
%     % this for-loop can be divided into three segments: norm
%     % computation of jth column and of jth with ith column, multiplying
%     % jth vector with this norm ratio and substraction.
%     % norm of jth column, and ith*jth column
%     for j = 1:i-1
%         [normjj,~] = innerprodTT(V{1,j},V{1,j});
%         [normij,~] = innerprodTT(V{1,i},V{1,j});
%         ratio = normij/normjj;
%         % multiply first core of jth column with ratio
%         V{1,j}{1,1} = V{1,j}{1,1}*ratio;
%         % subtract jth column from ith column, this can be done via
%         % addition of v(:,i) and -v(:,j)
%         V{1,i} = subTT(V{1,i},V{1,j});
%         % make function that outputs the ranks
%         rc(1) = 1;rc(d+1)=1;
%         for p = 1:d-1
%             c = size(V{1,i}{p,1});
%             rc(p+1) = c(3);
%         end
%         if max(rc) > 130
%             r = [1,2,2,1];
%             [V{1,i},rt] = roundTTr(V{1,i},r);
%         end
%     end    
%     Q1{1,i} = V{1,i};
%     [normQ,~] = innerprodTT(Q1{1,i},Q1{1,i});
%     Q1{1,i}{1,1} = Q1{1,i}{1,1}/sqrt(normQ);
% end
% 
% %% test Q1 for orthonormality
% [w,~] = innerprodTT(Q1{1,6},Q1{1,6});
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Gram-Schmidt algorithm for 1000 by 1000 matrix
% % Pieter van Klaveren
% % 01-04-2021
% clear all; clc; close all
% 
% A = rand(1000,1000)*5;
% a = size(A);
% v = zeros(size(A)); % orthogonal basis
% q = zeros(size(A)); % orthonormal basis
% for i = 1:a(2)
%     v(:,i) = A(:,i);
%     for j = 1:i-1
%         v(:,i) = v(:,i) - v(:,j)'*v(:,i)/(v(:,j)'*v(:,j))*v(:,j);
%     end    
%     q(:,i) = v(:,i)/norm(v(:,i));
% end
% 
% %% Gram-Schmidt algorithm in Tensor Train format for 20 by 20 matrix
% Row = [10,10,10];
% Col = [10,10,10];
% Aten = Mat2Tensor(A,Row,Col);
% eps = 0.01;
% ATT = Tensor2TT_SVD_eps(Aten,eps);
% ATTm = TT2TTm(ATT,Row,Col);
% clear Aten ATT
% 
% for i = 1:prod(Col)
%     % take ith vector out of TTm of A
%     % function to determine the code of indexing.
%     sizeATTm = size(ATTm);
%     e{1,1} = zeros(Col(1),1);
%     e{2,1} = zeros(Col(2),1);
%     e{3,1} = zeros(Col(3),1);
%     % determine e3
%     if i/Col(1)/Col(2) <= 1
%         e3 = 1;
%     elseif 1 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 2
%         e3 = 2;
%     elseif 2 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 3
%         e3 = 3;
%     elseif 3 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 4
%         e3 = 4;
%     elseif 4 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 5
%         e3 = 5;
%     elseif 5 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 6
%         e3 = 6;
%     elseif 6 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 7
%         e3 = 7;
%     elseif 7 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 8
%         e3 = 8;
%     elseif 8 < i/Col(1)/Col(2) && i/Col(1)/Col(2) <= 9
%         e3 = 9;
%     else 
%         e3 = 10;
%     end
%     e{3,1}(e3,1) = 1;
%     
%     % determine e2
%     if (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 1
%         e2 = 1;
%     elseif 1 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 2
%         e2 = 2;
%     elseif 2 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 3
%         e2 = 3;
%     elseif 3 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 4
%         e2 = 4;
%     elseif 4 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 5
%         e2 = 5;
%     elseif 5 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 6
%         e2 = 6;
%     elseif 6 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 7
%         e2 = 7;
%     elseif 7 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 8
%         e2 = 8;
%     elseif 8 < (i-(e3-1)*Col(1)*Col(2))/Col(1) && (i-(e3-1)*Col(1)*Col(2))/Col(1) <= 9
%         e2 = 9;
%     else 
%         e2 = 10;
%     end
%     e{2,1}(e2,1) = 1;
%     
%     % determine e1
%     if i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 1
%         e1 = 1;
%     elseif 1 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 2
%         e1 = 2;
%     elseif 2 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 3
%         e1 = 3;
%     elseif 3 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 4
%         e1 = 4;
%     elseif 4 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 5
%         e1 = 5;
%     elseif 5 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 6
%         e1 = 6;
%     elseif 6 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 7
%         e1 = 7;
%     elseif 7 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 8
%         e1 = 8;
%     elseif 8 < i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) && i-(e3-1)*Col(1)*Col(2)-(e2-1)*Col(1) <= 9
%         e1 = 9;
%     else 
%         e1 = 10;
%     end
%     e{1,1}(e1,1) = 1;
%     
%     % column pick from ATTm via e1,e2,e3
%     [d,~] = size(ATTm);
%     for k = 1:d
%         a = size(ATTm{k,1});
%         if k == d
%            a = [a,1]; 
%         end
%         V{1,i}{k,1} = reshape(permute(ATTm{k,1},[1,2,4,3]),[numel(ATTm{k,1})/a(3),a(3)])*e{k,1};
%         V{1,i}{k,1} = reshape(V{1,i}{k,1},[a(1),a(2),a(4)]);
%     end
%     
%     % this for-loop can be divided into three segments: norm
%     % computation of jth column and of jth with ith column, multiplying
%     % jth vector with this norm ratio and substraction.
%     % norm of jth column, and ith*jth column
%     for j = 1:i-1
%         [normjj,~] = innerprodTT(V{1,j},V{1,j});
%         [normij,~] = innerprodTT(V{1,i},V{1,j});
%         ratio = normij/normjj;
%         % multiply first core of jth column with ratio
%         V{1,j}{1,1} = V{1,j}{1,1}*ratio;
%         % subtract jth column from ith column, this can be done via
%         % addition of v(:,i) and -v(:,j)
%         V{1,i} = subTT(V{1,i},V{1,j});
%         % make function that outputs the ranks
%         rc(1) = 1;rc(d+1)=1;
%         for p = 1:d-1
%             c = size(V{1,i}{p,1});
%             rc(p+1) = c(3);
%         end
%         if max(rc) > 100
%             r = [1,8,8,1];
%             [V{1,i},rt] = roundTTr(V{1,i},r);
%         end
%     end    
%     Q1{1,i} = V{1,i};
%     [normQ,~] = innerprodTT(Q1{1,i},Q1{1,i});
%     Q1{1,i}{1,1} = Q1{1,i}{1,1}/sqrt(normQ);
% end
% 
% %% test Q1 for orthonormality
% [w,~] = innerprodTT(Q1{1,6},Q1{1,6});
end




