function TTm = ini_mat_eye(Mat1,Mat2,Row,Col)
% This function creates the covariance matrices via the kronecker product
% of two TT's having the combined dimension.

for i = 1:length(Row)
    if prod(Row(1:i)) == size(Mat2,1)
        break
    end
end

% Mat2: MxM
RowMat2 = Row(1:i);
ColMat2 = Col(1:i);
Mat2 = Mat2Tensor(Mat2,RowMat2,ColMat2);
eps = 1e-6;
Mat2 = Tensor2TT_SVD_eps(Mat2,eps);
Mat2 = TT2TTm(Mat2,RowMat2,ColMat2);
% Mat1: NxN
RowMat1 = Row(i+1:end);
ColMat1 = Col(i+1:end);
Mat1 = Mat2Tensor(Mat1,RowMat1,ColMat1);
eps = 1e-6;
Mat1 = Tensor2TT_SVD_eps(Mat1,eps);
Mat1 = TT2TTm(Mat1,RowMat1,ColMat1);
% TT-kronecker product
TTm.cores = [Mat2.cores;Mat1.cores];
TTm.dim = [Mat2.dim;Mat1.dim];

end