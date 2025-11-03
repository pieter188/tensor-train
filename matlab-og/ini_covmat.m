function Q = ini_covmat(data,Row,Col)
% function to compute initial process noise covariance matrix. The initial state covariance matrix is based on the 
% first frame

% INPUT: full video, downscaled
% OUTPUT: Initial state covariance matrix and initial process noise
% covariance matrix

%% from the data of the full video, but should be past video data
% if the amount of vertical pixels exceeds 40 than we select a random
% portion of the frame
if size(data,1) >= 40
    b = randi([1 size(data,1)-40],1);
    c = randi([1 size(data,2)-40],1);
    data1 = data(b:b+39,c:c+39,1,:);
else 
    data1 = data;
end
a = size(data1);
data1 = reshape(data1,[a(1)*a(2),a(4)]);
[W,~] = corrcoef(data1');
h = 1;
for j = 1:a(2)
    for i = 1:a(1)
        for jj = j:a(2)
            for ii = i:a(1)
                correlation(h) = W((jj-1)*a(1)+ii,(j-1)*a(1)+i);
                dist(h) = abs(i-ii) + abs(j-jj);
                h = h + 1;
            end
        end        
    end
end
for q = 1:a(1)+a(2)
    [~, j1] = find(dist==q-1);
    corr_mean(q) = mean(correlation(j1));
    dist_final(q) = q-1;
end

% compute bandwidth
k = find(corr_mean > 0.2);
bandwidth = max(k); %max(k)
%% create the Initial process covariance matrix
a = size(data);
Q1 = zeros(a(2)); %J/N
Q2 = zeros(a(1)); %I/M
for i = 1:bandwidth
    if i == 1
        Q1 = Q1 + eye(a(2));
        Q2 = Q2 + eye(a(1));
    else
        Q1 = Q1 + diag(ones(a(2)-i+1,1)*(1-i/(15+1)),-i+1)+diag(ones(a(2)-i+1,1)*(1-i/(15+1)),i-1);
        Q2 = Q2 + diag(ones(a(1)-i+1,1)*(1-i/(15+1)),-i+1)+diag(ones(a(1)-i+1,1)*(1-i/(15+1)),i-1);
    end
    
end

%% to TT-format
for i = 1:length(Row)
    if prod(Row(1:i)) == a(1)
        break
    end
end
% Q2: MxM / IxI
RowQ2 = Row(1:i);
ColQ2 = Col(1:i);
Q2 = Mat2Tensor(Q2,RowQ2,ColQ2);
eps = 1e-6;
Q2 = Tensor2TT_SVD_eps(Q2,eps);
Q2 = TT2TTm(Q2,RowQ2,ColQ2);

% Q1: NxN / JxJ
RowQ1 = Row(i+1:end);
ColQ1 = Col(i+1:end);
Q1 = Mat2Tensor(Q1,RowQ1,ColQ1);
eps = 1e-6;
Q1 = Tensor2TT_SVD_eps(Q1,eps);
Q1 = TT2TTm(Q1,RowQ1,ColQ1);



%% TT-kronecker product
Q.cores = [Q2.cores;Q1.cores];
Q.dim = [Q2.dim;Q1.dim];

end