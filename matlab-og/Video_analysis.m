%% Video analysis
% Pieter van Klaveren
% 25-02-2016
clear all; close all; clc

%% scaling of videos
% write a function that scales any resolution to any lower resolution
% video
% loading video
Video = VideoReader('original_video_TC.mp4');
% desired resolution
res = [480 720];
scaled_Video = scale_video(Video,res);
% display video on or off
playvideo = 0;
if playvideo == 1
    implay(scaled_Video,Video.FrameRate)
end

%% grayscale video
gray_Video = rgb2gray_video(scaled_Video);
% display video on or off
playvideo = 0;
if playvideo == 1
    implay(gray_Video,Video.FrameRate)
end
data = double(gray_Video);

%% Initial state covariance matrix
for i = 1:size(data,4)
    x{i} = reshape(data(:,:,1,i),[numel(data(:,:,1,i)),1]);
end
x_mean = (x{1}+x{2})/2;
% P0 = (x{1}-x_mean)*(x{1}-x_mean)';

sigma_p2 = zeros(1);
for i = 1:size(data,4)-1
    sigma_p2 = sigma_p2 + var(x{i+1}-x{i});
end
sigma_p2 = sigma_p2/size(data,4);
sigma_p = var(x{2}-x{1})/prod(res); %sqrt(sumsqr(x{2}-x{1})/prod(res));
% P0_s = sigma_p2*eye(prod(res));


s = var((x{2}-x{1}));

%% Initial process noise covariance matrix
% for i = 1:size(data,4)-1
%     PN_cov{i} = (x{i+1}-x{i})*(x{i+1}-x{i})';
% end
% PN = zeros(prod(res));
% for i = 1:size(data,4)-1
%     PN = PN + PN_cov{i};
% end
% PN = PN/(size(data,4));


% data = double(gray_Video);
% a = size(data);
% data = reshape(data,[a(1)*a(2),a(4)]);
% for k = 1:a(4)
%     w(:,k) = data(:,2)-data(:,1);
% end
% var(w(:))
%% Initial process noise covariance matrix
data = double(gray_Video);
data = data(100:200,300:400,1,:);
a = size(data);
data = reshape(data,[a(1)*a(2),a(4)]);
[W,P] = corrcoef(data');
% [W1,P1] = corr(data');
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


%% plot correlation vs pixels distance
figure
plot(dist_final(1:40),corr_mean(1:40),'k-','LineWidth',1.5)
grid on
ylim([0 1])
xlabel('Distance')
ylabel('Correlation')

%% create the Initial process covariance matrix
N = 720;
M = 480;
bandwidth = 15;
Q1 = eye(N);
Q2 = eye(M);
for i = 1:bandwidth
    Q1 = Q1 + diag(ones(N-i,1)*(1-i/(15+1)),-i)+diag(ones(N-i,1)*(1-i/(15+1)),i);
    Q2 = Q2 + diag(ones(M-i,1)*(1-i/(15+1)),-i)+diag(ones(M-i,1)*(1-i/(15+1)),i);
end

% to TT-format
% Q1: NxN
Row = [5 3 3 2 2 2 2];
Col = [5 3 3 2 2 2 2];
Q1 = Mat2Tensor(Q1,Row,Col);
eps = 1e-6;
Q1 = Tensor2TT_SVD_eps(Q1,eps);
Q1 = TT2TTm(Q1,Row,Col);

% Q1: MxM
Row = [5 3 2 2 2 2 2];
Col = [5 3 2 2 2 2 2];
Q2 = Mat2Tensor(Q2,Row,Col);
eps = 1e-6;
Q2 = Tensor2TT_SVD_eps(Q2,eps);
Q2 = TT2TTm(Q2,Row,Col);

% TT-kronecker product
Q.cores = [Q2.cores;Q1.cores];
Q.dim = [Q2.dim;Q1.dim];
