function Vf = combineTT(V,Cdim)
% This function combines multiple Tensor Trains into a Tensor Train matrix
% by concatenating the Tensor Trains.
% INPUT: V (stack of Tensor Train vector with the same row dimension),
% Cdim (Column dimension of the to be made Tensor Train matrix)
% OUTPUT: Tensor Train matrix, composed of the stack of Tensor Trains

d = size(V{1,1}.dim,1);
for i = 1:prod(Cdim)
    E = multi_index(Cdim,i);
    for k = 1:d
        a = V{1,i}.dim(k,:);
        b = size(E{k,1});
        Vm{1,i}.cores{k,1} = reshape(V{1,i}.cores{k,1},[a(1),a(2),1,a(3)]);
        Vm{1,i}.cores{k,1} = reshape(permute(Vm{1,i}.cores{k,1},[1,2,4,3]),...
            [a(1)*a(2)*a(3),1])*E{k,1}';
        Vm{1,i}.cores{k,1} = reshape(Vm{1,i}.cores{k,1},[a(1),a(2),a(3),b(1)]);
        Vm{1,i}.cores{k,1} = permute(Vm{1,i}.cores{k,1},[1,2,4,3]);
        Vm{1,i}.dim(k,:) = [a(1),a(2),b(1),a(3)];
    end
    [Vm{1,i},R,C] = TTm2TT(Vm{1,i});
end
% sum the vectors up in TT format
Vf = Vm{1,1};
for i = 2:prod(Cdim)
    Vf = addTT(Vf,Vm{1,i});
%     Vf = sitek(Vf,1);
    rank_threshold = 5;
    if max(Vf.dim(:,1)) > rank_threshold
        eps = 1e-4;
        [Vf,~] = roundTT2(Vf,eps); % kijk hier naar
    end
end
Vf = TT2TTm(Vf,R,C);
%% checked 04/05/2021: 11% error to original vectors, max 1e-6 innerprod on non-diagonals
end

% for i = 1:prod(Cdim)
%     E = multi_index(Cdim,i);
%     [d,~] = size(V{1,i});
%     for k = 1:d
%         a = size(V{1,i}{k,1});
%         if length(a) == 2
%            a = [a,1]; 
%         end
%         b = size(E{k,1});
%         Vm{1,i}{k,1} = reshape(V{1,i}{k,1},[a(1),a(2),1,a(3)]);
%         Vm{1,i}{k,1} = reshape(permute(Vm{1,i}{k,1},[1,2,4,3]),[a(1)*a(2)*a(3),1])*E{k,1}';
%         Vm{1,i}{k,1} = reshape(Vm{1,i}{k,1},[a(1),a(2),a(3),b(1)]);
%         Vm{1,i}{k,1} = permute(Vm{1,i}{k,1},[1,2,4,3]);
%     end
%     [Vm{1,i},R,C] = TTm2TT(Vm{1,i});
% end
% % sum the vectors up
% Vf = Vm{1,1};
% for i = 2:prod(Cdim)
%     Vf = addTT(Vf,Vm{1,i});
%     Vf = sitek(Vf,1);
%     rc(1) = 1; rc(d+1) = 1;
%     for p = 1:d-1
%         c = size(Vf{p,1});
%         rc(p+1) = c(3);
%     end
%     rank_threshold = 50;
%     if max(rc) > rank_threshold
%         eps = 1e-5;
%         [Vf,r] = roundTT(Vf,eps); % kijk hier naar
%     end
% end
% Vf = TT2TTm(Vf,R,C);