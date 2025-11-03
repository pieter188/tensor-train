function TTp = MatVecTT(TTm,TT)
% This function compute the TT of a matrix vector product given in TT(m)
% format.

d = size(TTm.dim,1);            % obtain the amount of cores d. 
for i = 1:d
    a = TTm.dim(i,:);           % a and b denote the size of the cores for
    b = TT.dim(i,:);            %  respectively ATT and xTT.  
    % FOR ATT: FROM: i1 x i2 x i3 x i4 TO: i1i2i4 x i3
    % FOR xTT: FROM: j1 x j2 x j3 TO: j2 x j1j3
    % BTT: i1i2i4 x j1j3
    TTp.cores{i,1} = reshape(permute(TTm.cores{i,1},[1,2,4,3]),[a(1)*a(2)*a(4),a(3)])*...
        reshape(permute(TT.cores{i,1},[2,1,3]),[b(2),b(1)*b(3)]);
    % FROM: i1i2i4 x j1j3 TO: i1 x i2 x i4 x j1 x j3
    TTp.cores{i,1} = reshape(TTp.cores{i,1},[a(1),a(2),a(4),b(1),b(3)]);
    % FROM: i1 x i2 x i4 x j1 x j3 TO: i1 x j1 x i2 x i4 x j3
    TTp.cores{i,1} = permute(TTp.cores{i,1},[1,4,2,3,5]);
    % FROM: i1 x j1 x i2 x i4 x j3 TO: i1j1 x i2 x i4j3
    TTp.cores{i,1} = reshape(TTp.cores{i,1},[a(1)*b(1),a(2),a(4)*b(3)]);
    TTp.dim(i,:) = [a(1)*b(1),a(2),a(4)*b(3)];
end
%% probably correct 04/05/2021
end



% [d,~] = size(TTm);             % obtain the amount of cores d. 
% for i = 1:d
%     a = size(TTm{i,1});        % a and b denote the size of the cores for
%     b = size(TT{i,1});        %  respectively ATT and xTT.
%     if length(a) == 3
%         a = [a,1];   % for the last core, the 1 isn't included via the size function so is done here manually
%     end
%     if length(b) == 2
%         b = [b,1];
%     end
%     
%     % FOR ATT: FROM: i1 x i2 x i3 x i4 TO: i1i2i4 x i3
%     % FOR xTT: FROM: j1 x j2 x j3 TO: j2 x j1j3
%     % BTT: i1i2i4 x j1j3
%     TTp{i,1} = reshape(permute(TTm{i,1},[1,2,4,3]),[a(1)*a(2)*a(4),a(3)])*...
%         reshape(permute(TT{i,1},[2,1,3]),[b(2),b(1)*b(3)]);
%     % FROM: i1i2i4 x j1j3 TO: i1 x i2 x i4 x j1 x j3
%     TTp{i,1} = reshape(TTp{i,1},[a(1),a(2),a(4),b(1),b(3)]);
%     % FROM: i1 x i2 x i4 x j1 x j3 TO: i1 x j1 x i2 x i4 x j3
%     TTp{i,1} = permute(TTp{i,1},[1,4,2,3,5]);
%     % FROM: i1 x j1 x i2 x i4 x j3 TO: i1j1 x i2 x i4j3
%     TTp{i,1} = reshape(TTp{i,1},[a(1)*b(1),a(2),a(4)*b(3)]);
%     
% end
% %% probably correct 30/10/2020

