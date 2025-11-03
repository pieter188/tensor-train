%% TNSRKF normal
% Pieter van Klaveren
% 25-02-2016
clear all; close all; clc
addpath('C:\Users\Pieter van Klaveren\Documents\Pieter van Klaveren\TU Delft\Master System and Control\SC52000-SC52045 Master Thesis\FinalThesis\simulation_data')

%% scaling of videos
% write a function that scales any resolution to any lower resolution
% video
% loading video
Video = VideoReader('original_video_TC.mp4');

% desired resolution
res = [9 16];

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

%% creating masked video
% fraction of missing pixels
alpha = 0.95;
% fraction of frames at the beginning that are not masked
beta = 0.01;
masked_Video = mask_video(gray_Video,alpha,beta);
% display video on or off
playvideo = 0;
if playvideo == 1
    implay(masked_Video,Video.FrameRate)
end

data = double(masked_Video);
a = size(data);

clear playvideo

%% define state space
n = numel(data(:,:,1,1));
A = eye(n);
% output data
y = reshape(data,[a(1)*a(2),a(4)]);

%% initial state and covariance matrices
% % state
% x_est_KF = zeros(a(1)*a(2),a(4));
% x_est_KF(:,1) = y(:,1);
% % covariance matrix
% P = eye(n);
% % process noise covariance matrix
% Q = eye(n);

% %% implement Kalman filter
% % this algorithm takes ages to run and S is singular which makes the kalman
% % algorithm unstable and causes the state estimate to diverge.
% for k = 1:(a(4)-1)
%     % computing C(k)
%     c = find(y(:,k+1)); 
%     C = zeros(length(c),n);
%     C(sub2ind(size(C),[1:length(c)]',c)) = 1;
%     % throwing away zeros in y(k)
%     y_atm = y(:,k+1);
%     y_atm = y_atm(y_atm~=0);
%     % time propagation
%     x_est_KF(:,k+1) = x_est_KF(:,k);
%     P = P + Q;
%     % measurement update
%     v = y_atm - C*x_est_KF(:,k+1);
%     S = C*P*C';
%     K = P*C'/S;
%     x_est_KF(:,k+1) = x_est_KF(:,k+1) + K*v;
%     P = P - K*C*P;
%     k
% end
% 
% %% implement partitioned kalman filter 
% % state
% x_est_PKF = zeros(a(1)*a(2),a(4));
% x_est_PKF(:,1) = y(:,1);
% % covariance matrix
% P = eye(n)*0.1;
% % process noise covariance matrix
% Q = eye(n)*5;
% for k = 1:(a(4)-1)
%     % computing C(k)
%     c = find(y(:,k+1)); 
%     C = zeros(length(c),n);
%     C(sub2ind(size(C),[1:length(c)]',c)) = 1;
%     % throwing away zeros in y(k)
%     y_atm = y(:,k+1);
%     y_atm = y_atm(y_atm~=0);
%     % time propagation
%     x_est_PKF(:,k+1) = x_est_PKF(:,k);
%     P = P + Q;
%     % measurement update
%     for j = 1:length(y_atm)
%         e = find(C(j,:));
%         v = y_atm(j) - x_est_PKF(e,k+1);
%         s = P(e,e);
%         K = P(:,e)/s;
%         x_est_PKF(:,k+1) = x_est_PKF(:,k+1) + K*v;
%         P = P - K*P(e,:);
%         j
%     end
%     k
% end
% 
% %% implement square-root kalman filter (with zeros in MU update)
% % initial square-root covariance matrix
% P = eye(n)*0.3+diag(ones(n-1,1)*0.1,-1)+diag(ones(n-2,1)*0.05,-2)+ ...
%     diag(ones(n-1,1)*0.1,1)+diag(ones(n-2,1)*0.05,2);
% L = chol(P); L = L';
% clear P
% % initial  square-root process noise covariance matrix
% Q = chol(eye(n)*5);
% % state
% x_SRKF = zeros(a(1)*a(2),a(4));
% x_SRKF(:,1) = y(:,1);
% 
% for k = 1:(a(4)-1)
%     % computing C(k)
%     c = find(y(:,k+1));
%     C = zeros(length(c),n);
%     C(sub2ind(size(C),[1:length(c)]',c)) = 1;
%     % throwing away zeros in y(k)
%     y_atm = y(:,k+1);
%     y_atm = y_atm(y_atm~=0);
%     if length(y_atm) == n
%         x_SRKF(:,k+1) = y_atm;
%     else
%         % time propagation
%         x_SRKF(:,k+1) = x_SRKF(:,k);
%         [~,L] = qr([L Q]');
%         L = L';
%         L = L(1:n,1:n);
%         % measurement update
%         v = y_atm - C*x_SRKF(:,k+1);
%         S = C*(L*L')*C';
%         K = (L*L')*C'/S;
%         x_SRKF(:,k+1) = x_SRKF(:,k+1) + K*v;
%         [~,L] = qr([zeros(length(c)), C*L; zeros(n,length(c)), L]');
%         L = L';
%         L = L((end-n+1):end,(end-n+1):end);
%     end
%     k
% end
% 
% % error analysis
% x_SRKF = reshape(x_SRKF,[res,301]);
% error = x_SRKF - reshape(double(gray_Video),[res,301]);
% error = permute(reshape(error,[prod(res),301]),[2,1]);
% ori = permute(reshape(gray_Video,[prod(res),301]),[2,1]);
% erSRKF = sqrt(sum(error.^2,2))./sqrt(sum(ori.^2,2));
% 
% 
% %% implement square-root kalman filter (without zeros in MU update)
% % initial square-root covariance matrix
% P = eye(n)*0.3+diag(ones(n-1,1)*0.1,-1)+diag(ones(n-2,1)*0.05,-2)+ ...
%     diag(ones(n-1,1)*0.1,1)+diag(ones(n-2,1)*0.05,2);
% L = chol(P); L = L';
% clear P
% % initial  square-root process noise covariance matrix
% Q = chol(eye(n)*5);
% % state
% x_SRKF2 = zeros(a(1)*a(2),a(4));
% x_SRKF2(:,1) = y(:,1);
% 
% for k = 1:(a(4)-1)
%     % computing C(k)
%     c = find(y(:,k+1));
%     C = zeros(length(c),n);
%     C(sub2ind(size(C),[1:length(c)]',c)) = 1;
%     % throwing away zeros in y(k)
%     y_atm = y(:,k+1);
%     y_atm = y_atm(y_atm~=0);
%     if length(y_atm) == n
%         x_SRKF2(:,k+1) = y_atm;
%     else
%         % time propagation
%         x_SRKF2(:,k+1) = x_SRKF2(:,k);
%         [~,L] = qr([L Q]');
%         L = L';
%         L = L(1:n,1:n);
%         % measurement update
%         v = y_atm - C*x_SRKF2(:,k+1);
%         S = C*(L*L')*C';
%         K = (L*L')*C'/S;
%         x_SRKF2(:,k+1) = x_SRKF2(:,k+1) + K*v;
%         [~,L] = qr([C*L; L]');
%         L = L';
%         L = L(length(c)+1:end,(end-n+length(c)+1):end);
%         L = [L,zeros(n,length(c))];
%     end
%     k
% end
% 
% % error analysis
% x_SRKF2 = reshape(x_SRKF2,[res,301]);
% error = x_SRKF2 - reshape(double(gray_Video),[res,301]);
% error = permute(reshape(error,[prod(res),301]),[2,1]);
% ori = permute(reshape(gray_Video,[prod(res),301]),[2,1]);
% erSRKF2 = sqrt(sum(error.^2,2))./sqrt(sum(ori.^2,2));
% 
% %% Implement partitioned square-root Kalman filter (with zeros in MU update)
% % initial square-root covariance matrix
% P = eye(n)*0.3+diag(ones(n-1,1)*0.1,-1)+diag(ones(n-2,1)*0.05,-2)+ ...
%     diag(ones(n-1,1)*0.1,1)+diag(ones(n-2,1)*0.05,2);
% L = chol(P); L = L';
% clear P
% % initial  square-root process noise covariance matrix
% Q = chol(eye(n)*5);
% % state
% x_PSRKF = zeros(a(1)*a(2),a(4));
% x_PSRKF(:,1) = y(:,1);
% for k = 1:(a(4)-1)
%     % computing C(k)
%     c = find(y(:,k+1));
%     C = zeros(length(c),n);
%     C(sub2ind(size(C),[1:length(c)]',c)) = 1;
%     % throwing away zeros in y(k)
%     y_atm = y(:,k+1);
%     y_atm = y_atm(y_atm~=0);
%     if length(y_atm) == n
%         x_PSRKF(:,k+1) = y_atm;
%     else
%         % time propagation
%         x_PSRKF(:,k+1) = x_PSRKF(:,k);
%         [~,L] = qr([L Q]');% A kan weg
%         L = L';
%         L = L(1:n,1:n);
%         % measurement update
%         for j = 1:length(y_atm)
%             e = find(C(j,:));
%             v = y_atm(j) - x_PSRKF(e,k+1);
%             s = L(e,:)*L(e,:)';
%             K = L*L(e,:)'/s;
%             x_PSRKF(:,k+1) = x_PSRKF(:,k+1) + K*v;
%             [~,L] = qr([zeros(1), L(e,:); zeros(n,1), L]');
%             L = L';
%             L = L((end-n+1):end,(end-n+1):end);
%             j
%         end
%     end
%     k
% end
% 
% % error analysis
% x_PSRKF = reshape(x_PSRKF,[res,301]);
% error = x_PSRKF - reshape(double(gray_Video),[res,301]);
% error = permute(reshape(error,[prod(res),301]),[2,1]);
% ori = permute(reshape(gray_Video,[prod(res),301]),[2,1]);
% erPSRKF = sqrt(sum(error.^2,2))./sqrt(sum(ori.^2,2));


%% Implement partitioned square-root Kalman filter (without zeros in MU update) via qr function
% initial square-root covariance matrix
L = eye(n);
% initial  square-root process noise covariance matrix
Q = eye(n);
% state
x = zeros(a(1)*a(2),a(4));
x(:,1) = y(:,1);
for k = 1:(a(4)-1)
    % computing C(k)
    c = find(y(:,k+1));
    C = zeros(length(c),n);
    C(sub2ind(size(C),[1:length(c)]',c)) = 1;
    % throwing away zeros in y(k)
    y_atm = y(:,k+1);
    y_atm = y_atm(y_atm~=0);
    if length(y_atm) == n
        x(:,k+1) = y_atm;
    else
        % time propagation
        x(:,k+1) = x(:,k);
        [~,L] = qr([L Q]');% A kan weg
        L = L';
        L = L(1:n,1:n);
        % measurement update
        for j = 1:length(y_atm)
            e = find(C(j,:));
            v = y_atm(j) - x(e,k+1);
            s = L(e,:)*L(e,:)';
            K = L*L(e,:)'/s;
            x(:,k+1) = x(:,k+1) + K*v;
            [~,L] = qr([L(e,:); L]');
            L = L';
            L = L(1+1:end,(end-n+1+1):end);
            L = [L,zeros(n,1)];
            j
        end
    end
    k
end

%% Implement partitioned square-root Kalman filter (without zeros in MU update) with own QR via MGS
% initial square-root covariance matrix
L = eye(n);
% initial  square-root process noise covariance matrix
Q = eye(n);
% state
xo = zeros(a(1)*a(2),a(4));
xo(:,1) = y(:,1);
for k = 1:(a(4)-1)
    % computing C(k)
    c = find(y(:,k+1));
    C = zeros(length(c),n);
    C(sub2ind(size(C),[1:length(c)]',c)) = 1;
    % throwing away zeros in y(k)
    y_atm = y(:,k+1);
    y_atm = y_atm(y_atm~=0);
    if length(y_atm) == n
        xo(:,k+1) = y_atm;
    else
        % time propagation
        xo(:,k+1) = xo(:,k);
        [~,L] = qr([L Q]');% A kan weg
        L = L';
        L = L(1:n,1:n);
        % measurement update
        for j = 1:length(y_atm)
            e = find(C(j,:));
            v = y_atm(j) - xo(e,k+1);
            s = L(e,:)*L(e,:)';
            K = L*L(e,:)'/s;
            xo(:,k+1) = xo(:,k+1) + K*v;
            v1 = [L(e,:); L]'; % orthogonal basis
            for i = 1:a(1)*a(2)
                v1(:,i) = P(:,i);
            end
            for i = 1:a(1)*a(2)
                q1(:,i) = v1(:,i)/sqrt(v1(:,i)'*v1(:,i));
                for j = i+1:a(1)*a(2)
                    v1(:,j) = v1(:,j) - (q1(:,i)'*v1(:,j))*q1(:,i);
                end
            end
            
            [~,L] = qr([L(e,:); L]');
            L = L';
            L = L(1+1:end,(end-n+1+1):end);
            L = [L,zeros(n,1)];
            j
        end
    end
    k
end

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

% % error analysis
% x_PSRKF2 = reshape(x_PSRKF2,[res,301]);
% error = x_PSRKF2 - reshape(double(gray_Video),[res,301]);
% error = permute(reshape(error,[prod(res),301]),[2,1]);
% ori = permute(reshape(gray_Video,[prod(res),301]),[2,1]);
% erPSRKF2 = sqrt(sum(error.^2,2))./sqrt(sum(ori.^2,2));
% 
% % error plot
% load('x_SRKF.mat')
% load('x_SRKF2.mat')
% load('x_PSRKF.mat')
% load('x_PSRKF2.mat')
% close all
% figure
% plot(erSRKF,'k')
% hold on
% plot(erSRKF2,'m')
% plot(erPSRKF,'b')
% plot(erPSRKF2,'r')
% hold off
% legend('SRKF','SRKF2','PSRKF','PSRKF2')
% xlabel('frame')
% ylabel('error')

% %% Analysis Square-root Kalman filter, 72p, 90 procent
% % video
% load('SRKF_xest2_72p_90.mat');
% SRKF_xest2_72p_90 = reshape(uint8(SRKF_xest2_72p_90),[72 128 1 301]);
% % implay(SRKF_xest_72p_90,30)
% % error per frame
% load('SRKF_xest2_72p_90.mat');
% load('SRKF_scaled_Video_72p_90.mat');
% SRKF_scaled_Video_72p_90 = double(SRKF_scaled_Video_72p_90);
% SRKF_xest2_72p_90 = reshape(SRKF_xest2_72p_90,[72 128 301]);
% SRKF_xest2_72p_90 = permute(SRKF_xest2_72p_90,[1 2 4 3]);
% error = SRKF_scaled_Video_72p_90 - SRKF_xest2_72p_90;
% for i=1:301
%     frame_error_SRKF2_90(i) = norm(error(:,:,1,i),'fro')/norm(SRKF_scaled_Video_72p_90(:,:,1,i));
% end
% clear i SRKF_scaled_Video_72p_90 error
% %% Analysis Square-root Kalman filter, 72p, 80 procent
% % video
% load('SRKF_xest_72p_80.mat');
% SRKF_xest_72p_80 = reshape(uint8(SRKF_xest_72p_80),[72 128 1 301]);
% % implay(SRKF_xest_72p_80,30)
% % error per frame
% load('SRKF_xest_72p_80.mat');
% load('SRKF_scaled_Video_72p_80.mat');
% SRKF_scaled_Video_72p_80 = double(SRKF_scaled_Video_72p_80);
% SRKF_xest_72p_80 = reshape(SRKF_xest_72p_80,[72 128 301]);
% SRKF_xest_72p_80 = permute(SRKF_xest_72p_80,[1 2 4 3]);
% error = SRKF_scaled_Video_72p_80 - SRKF_xest_72p_80;
% for i=1:301
%     frame_error_SRKF_80(i) = norm(error(:,:,1,i),'fro')/norm(SRKF_scaled_Video_72p_80(:,:,1,i));
% end
% clear i SRKF_scaled_Video_72p_80 error
% %% Analysis Kalman filter, 72p, 90 procent
% % video
% load('KF_xest_72p_90.mat');
% KF_xest_72p_90 = reshape(uint8(KF_xest_72p_90),[72 128 1 301]);
% % implay(KF_xest_72p_90,30)
% % error per frame
% load('KF_xest_72p_90.mat');
% load('KF_scaled_Video_72p_90.mat');
% KF_scaled_Video_72p_90 = double(KF_scaled_Video_72p_90);
% KF_xest_72p_90 = reshape(KF_xest_72p_90,[72 128 301]);
% KF_xest_72p_90 = permute(KF_xest_72p_90,[1 2 4 3]);
% error = KF_scaled_Video_72p_90 - KF_xest_72p_90;
% for i=1:301
%     frame_error_KF_90(i) = norm(error(:,:,1,i),'fro')/norm(KF_scaled_Video_72p_90(:,:,1,i));
% end
% clear i KF_scaled_Video_72p_90 error
% %% Analysis partitioned Kalman filter, 72p, 90 procent
% % video
% load('PKF_xest_72p_90.mat');
% PKF_xest_72p_90 = reshape(uint8(PKF_xest_72p_90),[72 128 1 301]);
% % implay(PKF_xest_72p_90,30)
% % error per frame
% load('PKF_xest_72p_90.mat');
% load('PKF_scaled_Video_72p_90.mat');
% PKF_scaled_Video_72p_90 = double(PKF_scaled_Video_72p_90);
% PKF_xest_72p_90 = reshape(PKF_xest_72p_90,[72 128 301]);
% PKF_xest_72p_90 = permute(PKF_xest_72p_90,[1 2 4 3]);
% error = PKF_scaled_Video_72p_90 - PKF_xest_72p_90;
% for i=1:301
%     frame_error_PKF_90(i) = norm(error(:,:,1,i),'fro')/norm(PKF_scaled_Video_72p_90(:,:,1,i));
% end
% clear i PKF_scaled_Video_72p_90 error
% %% Analysis partitioned square-root Kalman filter, 18p, 90 procent
% % video
% load('PSRKF_xest_18p_90.mat');
% PSRKF_xest_18p_90 = reshape(uint8(PSRKF_xest_18p_90),[18 32 1 301]);
% % implay(PKF_xest_72p_90,30)
% % error per frame
% load('PSRKF_xest_18p_90.mat');
% load('PSRKF_scaled_video_18p_90.mat');
% PSRKF_scaled_video_18p_90 = double(PSRKF_scaled_video_18p_90);
% PSRKF_xest_18p_90 = reshape(double(PSRKF_xest_18p_90),[18 32 301]);
% PSRKF_xest_18p_90 = permute(PSRKF_xest_18p_90,[1 2 4 3]);
% error = PSRKF_scaled_video_18p_90 - PSRKF_xest_18p_90;
% for i=1:301
%     frame_error_PSRKF_18p_90(i) = norm(error(:,:,1,i),'fro')/norm(PSRKF_scaled_video_18p_90(:,:,1,i));
% end
% clear i PSRKF_scaled_video_18p_90 error
% 
% %% plot
% close all
% figure
% plot(frame_error_KF_90,'k')
% hold on
% plot(frame_error_PKF_90,'b')
% % plot(frame_error_SRKF_90,'r--')
% plot(frame_error_SRKF2_90,'r')
% % plot(frame_error_SRKF_80,'m')
% % plot(frame_error_PSRKF_18p_90,'g')
% legend('KF90', 'PKF90', 'SRKF90')
% xlabel('frame')
% ylabel('error')

%% comparison of 
% load('x_0001')
% load('x_001')
% load('x_01')
% load('x_1')
% load('x_10')
% load('x_100')
% load('x_1000')
% gray_Video = reshape(double(gray_Video),[res,301]);
% x_0001 = reshape(x_0001,[res,301]);
% x_001 = reshape(x_001,[res,301]);
% x_01 = reshape(x_01,[res,301]);
% x_1 = reshape(x_1,[res,301]);
% x_10 = reshape(x_10,[res,301]);
% x_100 = reshape(x_100,[res,301]);
% x_1000 = reshape(x_1000,[res,301]);
% 
% error0001 = gray_Video - x_0001;
% error001 = gray_Video - x_001;
% error01 = gray_Video - x_01;
% error1 = gray_Video - x_1;
% error10 = gray_Video - x_10;
% error100 = gray_Video - x_100;
% error1000 = gray_Video - x_1000;
% 
% frame_error0001 = zeros(1,301);
% frame_error001 = zeros(1,301);
% frame_error01 = zeros(1,301);
% frame_error1 = zeros(1,301);
% frame_error10 = zeros(1,301);
% frame_error100 = zeros(1,301);
% frame_error1000 = zeros(1,301);
% for i=1:301
%     frame_error0001(i) = norm(error0001(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
%     frame_error001(i) = norm(error001(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
%     frame_error01(i) = norm(error01(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
%     frame_error1(i) = norm(error1(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
%     frame_error10(i) = norm(error10(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
%     frame_error100(i) = norm(error100(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
%     frame_error1000(i) = norm(error1000(:,:,i),'fro')/norm(gray_Video(:,:,i),'fro');
% end
% close all
% figure
% hold on
% plot(frame_error0001,'k')
% plot(frame_error001,'m')
% plot(frame_error01,'g')
% plot(frame_error1,'b')
% plot(frame_error10,'r')
% plot(frame_error100,'y')
% plot(frame_error1000,'k--')
% hold off
% legend('x0001','x001','x01','x1','x10','x100','x1000')
% xlabel('frame')
% ylabel('error')
