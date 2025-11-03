%% main document
% Pieter van Klaveren
% 25-02-2016
clear all; close all; clc

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
data = double(masked_Video(:,:,:,:));
a = size(data);

%% Preparation tensor networked square-root Kalman filter
% define state space in TT-format
n = numel(data(:,:,1,1));
% partitioning of row and column of A
Row = [5 3 2 2 2 2 2 5 3 3 2 2 2 2];
Col = [5 3 2 2 2 2 2 5 3 3 2 2 2 2];
A1 = eye(size(data,2));
A2 = eye(size(data,1));
ATTm = ini_mat_eye(A1,A2,Row,Col);

% C matrix
C1 = eye(size(data,2));
C2 = eye(size(data,1));
CTTm = ini_mat_eye(C1,C2,Row,Col);

% initial square-root covariance matrix
L1 = eye(size(data,2));
L2 = eye(size(data,1));
LTTm = ini_mat_eye(L1,L2,Row,Col);

% process noise covariance matrix
QTTm = ini_covmat(data,Row,Col);

% output data
y = reshape(data,[a(1)*a(2),a(4)]);
% state
x_est = zeros(a(1)*a(2),a(4));
x_est(:,1) = y(:,1);
xten = Vec2Tensor(x_est(:,1),Row);
xTT{1,1} = Tensor2TT_SVD_eps(xten,eps);
clear xten xest

%% Implementation tensor networked square-root Kalman filter
for k = 1:a(4)-1
    k
    % check whether y_atm is sparse
    [ii,~,~] = find(y(:,k+1));
    if size(ii,1) == size(y(:,k+1),1)
        yten = Vec2Tensor(y(:,k+1),Row);
        eps = 1e-6;
        yTT = Tensor2TT_SVD_eps(yten,eps);
        xTT{1,k+1} = yTT;
        xTT{1,k+1} = sitek(xTT{1,k+1},1);
    else
        %% time propagation
        xTT{1,k+1} = xTT{1,k};
        LTTm = concatTTm(LTTm,QTTm);
        % take the transpose because we are computing QR
        LTTm = TransposeTTm(LTTm);
        % bring LTTm into site-1 canonical form
        [LTT,RowL,ColL] = TTm2TT(LTTm);        
        LTTc = sitek(LTT,1);
        LTTm = TT2TTm(LTTc,RowL,ColL);
        % QR decomposition of LTTm
        % this for-loop extracts the column vectors from LTTm
        for i = 1:prod(Col)
            % take ith vector out of TTm of A
            % function to determine the code of indexing.
            e = multi_index(Col,i);
            % column pick from ATTm via e1,e2,e3
            d = size(LTTm.dim,1);
            for m = 1:d
                V.cores{m,1} = reshape(permute(LTTm.cores{m,1},[1,2,4,3]),...
                    [numel(LTTm.cores{m,1})/LTTm.dim(m,3),LTTm.dim(m,3)])*e{m,1};
                V.cores{m,1} = reshape(V.cores{m,1},...
                    [LTTm.dim(m,1),LTTm.dim(m,2),LTTm.dim(m,4)]);
                V.dim(m,:) = [LTTm.dim(m,1),LTTm.dim(m,2),LTTm.dim(m,4)];
            end
            V_tp{1,i} = V;
        end   
        % NOTE: ALTHOUGH LTTm IS IN SITE-1 CANONICAL FORM, THE VECTOR
        % SUBSTRACTED FROM LTTm ARE NOT ANY LONGER IN SITE-1 CANONICAL FORM.

        % this for-loop can be divided into three segments: norm
        % computation of jth column and of jth with ith column, multiplying
        % jth vector with this norm ratio and substraction.
        % norm of jth column, and ith*jth column
        for i = 1:prod(Col)
            % normalize Q via site-1 orthogonalization and division by first
            % core norm
            Q_tp{1,i} = V_tp{1,i};
            Q_tp{1,i} = sitek(Q_tp{1,i},1);
            normQ1 = norm(reshape(Q_tp{1,i}.cores{1,1},[numel(Q_tp{1,i}.cores{1,1}),1]),'fro');
            Q_tp{1,i}.cores{1,1} = Q_tp{1,i}.cores{1,1}/normQ1;
            for j = i+1:prod(Col)
                normij = innerprodTT(Q_tp{1,i},V_tp{1,j});
                QQ = Q_tp{1,i};
                QQ.cores{1,1} = normij*QQ.cores{1,1};
                V_tp{1,j} = subTT(V_tp{1,j},QQ);
                V_tp{1,j} = sitek(V_tp{1,j},1);
                % check rank and truncate
                rank_threshold = 10;
                if max(V_tp{1,j}.dim(1,:)) > rank_threshold
                    eps = 1*10^(-5);
%                     [V_tp{1,j},~] = roundTT(V_tp{1,j},eps); % kijk hier naar
                    [V_tp{1,j},~] = roundTT2(V_tp{1,j},eps); % kijk hier naar
                end
            end    
        end
        % vectors in Q1 added together into one TTm
        Q = combineTT(Q_tp,Col);
        % it is checked that after this function Q^T*Q is still identity,
        % an error up to e-6 is observed.
        % compute R
        R = MatMatTT(TransposeTTm(Q),LTTm); 
        [R,RowR,ColR] = TTm2TT(R);
        eps = 1e-6;
        [R,~] = roundTT2(R,eps);
        LTTm = TT2TTm(R,RowR,ColR);
        LTTm = TransposeTTm(LTTm);
        % Time propagation has been checked and is good. The new LTTm has
        % error of e-6 to the original R1.
        % This new LTTm has the (original) desired dimensions equal to
        % prod(res)xprod(res).
        % ----------------------------------------------------------------
        % This section is to check whether the measurement update is
        % working correctly.
%         L_tp = TTm2TT(LTTm);
%         L_tp = TT2Tensor(L_tp);
%         L_tp = Tensor2Mat(L_tp,RowR,ColR);
        
        %% Partitioned measurement update
        [ii,~,y_sp1] = find(y(:,k+1));
        % In this massive for-loop x{k+1} and LTTm are being updated
        % according to each single measurement obtained. This means that
        % the new LTTm must remean the same dimension to avoid unnecessary
        % growth and that each new update uses LTTm from the previous
        % update.
        % ------------------------------------------------------------
        % This section takes the correct vector out of CTTm to compute
        % residual vector v.
        % C matrix
        C1 = eye(size(data,2));
        C2 = eye(size(data,1));
        CTTm = ini_mat_eye(C1,C2,Row,Col);
        for l = 1:length(ii)
            e = multi_index(Col,ii(l));
            % column pick from CTTm via e1,e2,e3
            d = size(CTTm.dim,1);
            for m = 1:d
                a = CTTm.dim(m,:);
                CTT{1,l}.cores{m,1} = reshape(permute(CTTm.cores{m,1},[1,2,4,3]),...
                    [numel(CTTm.cores{m,1})/a(3),a(3)])*e{m,1};
                CTT{1,l}.cores{m,1} = reshape(CTT{1,l}.cores{m,1},[a(1),a(2),a(4)]);
                CTT{1,l}.dim(m,:) = [a(1),a(2),a(4)];
            end
        end
        CTTm = combineTT(CTT,[ones(1,length(Row)-1),length(ii)]);
        CTTm = TransposeTTm(CTTm);
        [CTTm,RowC,ColC] = TTm2TT(CTTm);
        eps = 1e-5;
        CTTm = roundTT2(CTTm,eps);
        CTTm = TT2TTm(CTTm,RowC,ColC);
        if length(y_sp1) == 1
            for m = 1:d
                y_sp.cores{m,1} = 1; 
                y_sp.dim(m,:) = [1 1 1];
                if m == d
                    y_sp.cores{m,1} = y_sp1;
                end
            end
        else
            y_sp = y_sp1;
            y_sp = Vec2Tensor(y_sp,[ones(1,length(Row)-1),length(ii)]);
            eps = 1e-5;
            y_sp = Tensor2TT_SVD_eps(y_sp,eps);
        end
        
        v = subTT(y_sp,MatVecTT(CTTm,xTT{1,k+1}));
        eps = 1e-5;
        v = roundTT2(v,eps);
        % ------------------------------------------------------------
        % This section computes the residual covariance matrix, the
        % kalman gain and the (measurement) updated x{k+1}
        s = MatMatTT(TransposeTTm(LTTm),TransposeTTm(CTTm));
        [s,Rows,Cols] = TTm2TT(s);
        s = roundTT2(s,eps);
        s = TT2TTm(s,Rows,Cols);
        s = MatMatTT(CTTm,MatMatTT(LTTm,s));
        s = TTm2TT(s);
        s = roundTTr(s,[1,1,1,1,1,1]);
        s = TT2Tensor(s);
        if length(y_sp1) == 1
            s = s;
        else
            s = Tensor2Mat(s,[ones(1,length(Row)-1),length(ii)],...
                [ones(1,length(Row)-1),length(ii)]);
        end
        KTT = MatMatTT(LTTm,MatMatTT(TransposeTTm(LTTm),TransposeTTm(CTTm)));
        KTT = MatVecTT(KTT,v);
        KTT = roundTT2(KTT,eps);
        KTT.cores{1,1} = KTT.cores{1,1}/max(max(s));       
        xTT{1,k+1} = addTT(xTT{1,k+1},KTT);
        xTT{1,k+1} = sitek(xTT{1,k+1},1);
        % ------------------------------------------------------------
        % This section takes the vectors out of LTTm and stores them
        % into L_vec. In total there are prod(res) vectors, since this
        % is the size of LTTm. Save the vectors for later.
        LTTm = TransposeTTm(LTTm);
        for i = 1:prod(Col)
            % take ith vector out of TTm of A
            % function to determine the code of indexing.
            e = multi_index(Col,i);
            % column pick from ATTm via e1,e2,e3
            d = size(LTTm.dim,1);
            for m = 1:d
                a = LTTm.dim(m,:);
                V_mu{1,i+length(ii)}.cores{m,1} = reshape(permute(LTTm.cores{m,1},[1,2,4,3]),...
                    [numel(LTTm.cores{m,1})/a(3),a(3)])*e{m,1};
                V_mu{1,i+length(ii)}.cores{m,1} = reshape(V_mu{1,i+length(ii)}.cores{m,1},...
                    [a(1),a(2),a(4)]);
                V_mu{1,i+length(ii)}.dim(m,:) = [a(1),a(2),a(4)];
            end
        end
        L_vec = V_mu;
        D = MatMatTT(TransposeTTm(LTTm),TransposeTTm(CTTm));
        for i = 1:length(ii)
            % take ith vector out of TTm of A
            % function to determine the code of indexing.
            e = multi_index([ones(1,length(Row)-1),length(ii)],i);
            % column pick from ATTm via e1,e2,e3
            d = size(LTTm.dim,1);
            for m = 1:d
                a = D.dim(m,:); % size(LTTm{m,1})
                L_vec{1,i}.cores{m,1} = reshape(permute(D.cores{m,1},[1,2,4,3]),...
                    [numel(D.cores{m,1})/a(3),a(3)])*e{m,1};
                L_vec{1,i}.cores{m,1} = reshape(L_vec{1,i}.cores{m,1},...
                    [a(1),a(2),a(4)]);
                L_vec{1,i}.dim(m,:) = [a(1),a(2),a(4)];
            end
        end
        V_mu = L_vec;
        % ------------------------------------------------------------
        % This section computes Q via the Modified Gram-Schmidt
        % orthogonalization by using the vectors stored in V and making
        % them orthogonal to (C*L)'. Because the QR decomposition is
        % taken of Lnew' = [C*L;L]'. Say that n = prod(res). Then Lnew
        % is nxn+1, Q will then be nxn and R will be nxn+1. Thus only
        % the first n-1 vectors in V will be used to compute Q.         
        for i = 1:prod(Col)
            Q_mu{1,i} = V_mu{1,i};                    
            % normalize Q via site-1 orthogonalization and division by first
            % core norm
            Q_mu{1,i} = sitek(Q_mu{1,i},1);
            normQ1 = norm(reshape(Q_mu{1,i}.cores{1,1},[numel(Q_mu{1,i}.cores{1,1}),1]),'fro');
            Q_mu{1,i}.cores{1,1} = Q_mu{1,i}.cores{1,1}/normQ1;
            for j = i+1:prod(Col)
                normij = innerprodTT(Q_mu{1,i},V_mu{1,j});   
                QQ = Q_mu{1,i};
                QQ.cores{1,1} = normij*QQ.cores{1,1};
                V_mu{1,j} = subTT(V_mu{1,j},QQ);
                V_mu{1,j} = sitek(V_mu{1,j},1);
                % check rank and truncate
                rank_threshold = 10;
                if max(V_mu{1,j}.dim(1,:)) > rank_threshold
                    eps = 1*10^(-5);
                    [V_mu{1,j},~] = roundTT2(V_mu{1,j},eps);
                end
            end    
        end
        % Combine n orthonormal vectors into one matrix.
        Q = combineTT(Q_mu,Col);
        [Q,RowQ,ColQ] = TTm2TT(Q);
%             Q = sitek(Q,1);
        Q = roundTT2(Q,eps);
        Q = TT2TTm(Q,RowQ,ColQ);
        % ------------------------------------------------------------
        % This section creates LTTm = [C*L;L;O]. This is checked.
        for i = 1:d
            zz.cores{i,1} = zeros(1,Row(i),1);
            zz.dim(i,:) = [1,Row(i),1];
        end
        for i = prod(Col)+length(ii):prod(Col)*2
            L_vec{1,i} = zz;
        end
        Colz = Col;
        Colz(end) = Colz(end)*2;
        LTTm = combineTT(L_vec,Colz);
        [LTTm,ColL,RowL] = TTm2TT(LTTm);
        LTTm = roundTT2(LTTm,eps);
        LTTm = TT2TTm(LTTm,ColL,RowL);
        % ------------------------------------------------------------          
        % This section computes L (= R^T)
        R = MatMatTT(TransposeTTm(Q),LTTm);
        [R,RowR,ColR] = TTm2TT(R);
%             R = sitek(R,1);
        eps = 1e-4;
        [R,~] = roundTT2(R,eps);
        eps = 1e-5;
        R = TT2TTm(R,RowR,ColR);
        % ------------------------------------------------------------
        % This section takes R2 out of R by taking the 2-n+1 column of
        % R
        for i = length(ii)+1:prod(RowR)+length(ii)
            % function to determine the code of indexing.
            e = multi_index(ColR,i);
            % column pick from ATTm via e1,e2,e3
            d = size(R.dim,1);
            for m = 1:d
                a = R.dim(m,:);%size(R{m,1});
                R_col{1,i-length(ii)}.cores{m,1} = reshape(permute(R.cores{m,1},...
                    [1,2,4,3]),[numel(R.cores{m,1})/a(3),a(3)])*e{m,1};
                R_col{1,i-length(ii)}.cores{m,1} = reshape(R_col{1,i-length(ii)}.cores{m,1},...
                    [a(1),a(2),a(4)]);
                R_col{1,i-length(ii)}.dim(m,:) = [a(1),a(2),a(4)];
            end
        end
        R2 = combineTT(R_col,ColL);
        [R2,RowR2,ColR2] = TTm2TT(R2);
%             R2 = sitek(R2,1);
        [R2,~] = roundTT2(R2,eps);
        R2 = TT2TTm(R2,RowR2,ColR2);
        R2 = TransposeTTm(R2);
        % ------------------------------------------------------------
        % Add a zero column to the end of R2 and remove the first
        % column -> R22
        for i = length(ii)+1:prod(RowR2)
            % function to determine the code of indexing.
            e = multi_index(ColR2,i);
            % column pick from ATTm via e1,e2,e3
            d = size(R2.dim,1);
            for m = 1:d
                a = R2.dim(m,:);
                R2_col{1,i-length(ii)}.cores{m,1} = reshape(permute(R2.cores{m,1},...
                    [1,2,4,3]),[numel(R2.cores{m,1})/a(3),a(3)])*e{m,1};
                R2_col{1,i-length(ii)}.cores{m,1} = reshape(R2_col{1,i-length(ii)}.cores{m,1},...
                    [a(1),a(2),a(4)]);
                R2_col{1,i-length(ii)}.dim(m,:) = [a(1),a(2),a(4)];
            end
        end
        for i = prod(RowR2)-length(ii)+1:prod(RowR2)
            R2_col{1,i} = zz;
        end
        LTTm = combineTT(R2_col,ColR2);
        [LTTm,RowL,ColL] = TTm2TT(LTTm);
        eps = 1e-5;
        [LTTm,~] = roundTT2(LTTm,eps);
        LTTm = TT2TTm(LTTm,RowL,ColL);
%         LLL{1,k} = LTTm;
%             L22_test = TTm2TT(LTTm);
%             L22_test = TT2Tensor(L22_test);
%             L22_test = Tensor2Mat(L22_test,[3,4,4,3],[3,4,4,3]);
%             [QQQ,RRR] = qr([LC,L_tp]);
%             RRR = RRR'; RRR = [RRR(2:end,2:end),zeros(144,1)];
    end
end