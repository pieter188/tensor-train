function C = addTT(A,B)
% This function computes the addition of 2 tensor train. This is
% done according to: Tensor Train Decompositions by Oseledets, 2011.
% input: A and B
% output: C where C = A + B

d = size(A.dim,1);
for i = 1:d
    % reshape cores into a matrix such that the right rank is with the
    % index
    A.cores{i,1} = reshape(A.cores{i,1},[A.dim(i,1),A.dim(i,2)*A.dim(i,3)]);
    B.cores{i,1} = reshape(B.cores{i,1},[B.dim(i,1),B.dim(i,2)*B.dim(i,3)]);

    % compose C
    if i ==1
        C.cores{i,1} = [A.cores{i,1}, B.cores{i,1}];
        C.cores{i,1} = reshape(C.cores{i,1},[A.dim(i,1),A.dim(i,2),...
            A.dim(i,3)+B.dim(i,3)]);
        C.dim(i,:) = [A.dim(i,1),A.dim(i,2),A.dim(i,3)+B.dim(i,3)];
    elseif i == d
        C.cores{i,1} = [A.cores{i,1}; B.cores{i,1}];
        C.cores{i,1} = reshape(C.cores{i,1},[A.dim(i,1)+B.dim(i,1),...
            A.dim(i,2),A.dim(i,3)]);
        C.dim(i,:) = [A.dim(i,1)+B.dim(i,1),A.dim(i,2),A.dim(i,3)];
    else
        C.cores{i,1} = [A.cores{i,1}, zeros(A.dim(i,1),B.dim(i,2)*B.dim(i,3));
                    zeros(B.dim(i,1),A.dim(i,2)*A.dim(i,3)), B.cores{i,1}];
        C.cores{i,1} = reshape(C.cores{i,1},[A.dim(i,1)+B.dim(i,1),...
            A.dim(i,2),A.dim(i,3)+B.dim(i,3)]);
        C.dim(i,:) = [A.dim(i,1)+B.dim(i,1),A.dim(i,2),A.dim(i,3)+B.dim(i,3)];
    end
end
%% checked 01/05/2021
end

% [d,~] = size(A);
% for i = 1:d
%     % obtain size of each single core
%     a = size(A{i,1});
%     b = size(B{i,1});
%     if length(a) == 2
%        a = [a,1];  
%     end
%     if length(b) == 2
%         b = [b,1];
%     end
%     % reshape cores into a matrix such that the right rank is with the
%     % index
%     A{i,1} = reshape(A{i,1},[a(1),a(2)*a(3)]);
%     B{i,1} = reshape(B{i,1},[b(1),b(2)*b(3)]);
% 
%     % compose C
%     if i ==1
%         C{i,1} = [A{i,1}, B{i,1}];
%         C{i,1} = reshape(C{i,1},[a(1),a(2),a(3)+b(3)]);
%     elseif i == d
%         C{i,1} = [A{i,1}; B{i,1}];
%         C{i,1} = reshape(C{i,1},[a(1)+b(1),a(2),a(3)]);
%     else
%         C{i,1} = [A{i,1}, zeros(a(1),b(2)*b(3));
%                     zeros(b(1),a(2)*a(3)), B{i,1}];
%         C{i,1} = reshape(C{i,1},[a(1)+b(1),a(2),a(3)+b(3)]);
%     end
% end