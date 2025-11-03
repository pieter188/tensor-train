function TTp = MatMatTT(TTm1,TTm2)
% This function compute the TT of a matrix vector product given in TT(m)
% format.
% INPUT: Tensor-Train matrix 1, Tensor Train matrix 2
% OUTPUT: Product of TTm1 and TTm2

d = size(TTm1.dim,1);             % obtain the amount of cores d. 
for i = 1:d
    a = TTm1.dim(i,:);          % a and b denote the size of the cores for
    b = TTm2.dim(i,:);          %  respectively ATT and xTT.
    % FOR ATT: FROM: i1 x i2 x i3 x i4 TO: i1i2i4 x i3
    % FOR xTT: FROM: j1 x j2 x j3 x j4 TO: j2 x j1j3j4
    % BTT: i1i2i4 x j1j3j4
    TTp.cores{i,1} = reshape(permute(TTm1.cores{i,1},[1,2,4,3]),[a(1)*a(2)*a(4),a(3)])*...
        reshape(permute(TTm2.cores{i,1},[2,1,3,4]),[b(2),b(1)*b(3)*b(4)]);
    % FROM: i1i2i4 x j1j3j4 TO: i1 x i2 x i4 x j1 x j3 x j4
    TTp.cores{i,1} = reshape(TTp.cores{i,1},[a(1),a(2),a(4),b(1),b(3),b(4)]);
    % FROM: i1 x i2 x i4 x j1 x j3 TO: i1 x j1 x i2 x i4 x j3
    TTp.cores{i,1} = permute(TTp.cores{i,1},[1,4,2,5,3,6]);
    % FROM: i1 x j1 x i2 x i4 x j3 TO: i1j1 x i2 x i4j3
    TTp.cores{i,1} = reshape(TTp.cores{i,1},[a(1)*b(1),a(2),b(3),a(4)*b(4)]);
    TTp.dim(i,:) = [a(1)*b(1),a(2),b(3),a(4)*b(4)];
end
%% checked 04/05/2021
end

% 
% [d,~] = size(TTm1);             % obtain the amount of cores d. 
% for i = 1:d
%     a = size(TTm1{i,1});        % a and b denote the size of the cores for
%     b = size(TTm2{i,1});        %  respectively ATT and xTT.
%     if length(a) == 3
%         a = [a,1];   % for the last core, the 1 isn't included via the size function so is done here manually
%     end
%     if length(b) == 3
%         b = [b,1];
%     end
%     % FOR ATT: FROM: i1 x i2 x i3 x i4 TO: i1i2i4 x i3
%     % FOR xTT: FROM: j1 x j2 x j3 x j4 TO: j2 x j1j3j4
%     % BTT: i1i2i4 x j1j3j4
%     TTp{i,1} = reshape(permute(TTm1{i,1},[1,2,4,3]),[a(1)*a(2)*a(4),a(3)])*...
%         reshape(permute(TTm2{i,1},[2,1,3,4]),[b(2),b(1)*b(3)*b(4)]);
%     % FROM: i1i2i4 x j1j3j4 TO: i1 x i2 x i4 x j1 x j3 x j4
%     TTp{i,1} = reshape(TTp{i,1},[a(1),a(2),a(4),b(1),b(3),b(4)]);
%     % FROM: i1 x i2 x i4 x j1 x j3 TO: i1 x j1 x i2 x i4 x j3
%     TTp{i,1} = permute(TTp{i,1},[1,4,2,5,3,6]);
%     % FROM: i1 x j1 x i2 x i4 x j3 TO: i1j1 x i2 x i4j3
%     TTp{i,1} = reshape(TTp{i,1},[a(1)*b(1),a(2),b(3),a(4)*b(4)]);
%     
% end
% %% Checked on 26/03/2021


%     % FOR ATT: FROM: i1 x i2 x i3 x i4 TO: i1i2i4 x i3
%     % FOR xTT: FROM: j1 x j2 x j3 TO: j2 x j1j3
%     BTT{i,1} = reshape(permute(ATT{i,1},[4,2,1,3]),[a(4)*a(2)*a(1),a(3)])*...
%         reshape(permute(xTT{i,1},[2,3,1]),[b(2),b(3)*b(1)]);
%     B2 = reshape(permute(ATT{i,1},[4,2,1,3]),[a(1)*a(2)*a(4),a(3)])*...
%         reshape(permute(xTT{i,1},[2,3,1]),[b(2),b(1)*b(3)]);
%     % FROM: i1i2i3 x j1j3 TO: j3i1i2i4 x j1
%     BTT{i,1} = reshape(BTT{i,1},[a(4)*a(2)*a(1)*b(3),b(1)]);
%     % FROM: j3i1i2i4 x j1 TO: j1 x j3i1i2i4
%     BTT{i,1} = permute(BTT{i,1},[2,1]);
%     % FROM: j1 x j3i1i2i4 TO: i4j1 x j3i1i2
%     BTT{i,1} = reshape(BTT{i,1},[a(4)*b(1),a(1)*a(2)*b(3)]);
%     % FROM: i4j1 x j3i1i2 TO: j3i1i2 x i4j1
%     BTT{i,1} = permute(BTT{i,1},[2,1]);
%     % FROM: j3i1i2 x i4j1 TO: j3i1i2 x j1 x i4
%     BTT{i,1} = reshape(BTT{i,1},[a(1)*a(2)*b(3),b(1),a(4)]);
%     % FROM: j3i1i2 x j1 x i4 TO: j1 x j3i1i2 x i4
%     BTT{i,1} = permute(BTT{i,1},[2,1,3]);
%     % FROM: j1 x j3i1i2 x i4 TO: j1 x i2 x j3i1 x i4
%     BTT{i,1} = reshape(BTT{i,1},[b(1),a(2),a(1)*b(3),a(4)]);
%     % FROM: j1 x i2 x j3i1 x i4 TO: j1 x i2 x i1 x j3 x i4
%     BTT{i,1} = reshape(BTT{i,1},[b(1),a(2),a(1),b(3),a(4)]);
%     % FROM: j1 x i2 x i1 x j3 x i4 TO: j1 x i1 x i2 x j3 x i4
%     BTT{i,1} = permute(BTT{i,1},[1,3,2,4,5]);
%     % FROM: j1 x i1 x i2 x j3 x i4 TO: i1j1 x i2 x j3 x i4
%     BTT{i,1} = reshape(BTT{i,1},[b(1)*a(1),a(2),b(3),a(4)]);
%     % FROM: i1j1 x i2 x j3 x i4 TO: i1j1 x i2 x i4j3
%     BTT{i,1} = reshape(BTT{i,1},[b(1)*a(1),a(2),b(3)*a(4)]);

