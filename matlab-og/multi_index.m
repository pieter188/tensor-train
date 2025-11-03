function E = multi_index(I,i)
% input: Index I, current combined value i
% output: multi-index ordening: [i1, i2, ..., id]

num = i;
for j = length(I):-1:1
    for s = 1:I(j)
        if s-1 < num/prod(I(1:j-1)) && num/prod(I(1:j-1)) <= s
            e(j) = s;
            num = num - (e(j)-1)*prod(I(1:j-1));
        end
    end
    E{j,1} = zeros(I(j),1);
    E{j,1}(e(j)) = 1;
end
%% checked 03/05/2021
end