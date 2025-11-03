function Row = row_part(n)
% input: the amount of indices
% output: a vector containing the partitioning
% ik wil een functie die aan de hand van de grootte van een dimensie, deze
% dimensie opsplitst in verschillende delen. De splitsingen kunnen
% bijvoorbeeld tussen 2 en 5 liggen.

X0 = [3;3;3;3;3;3;3;3];
Aeq = [1,1,1,1,1,1,1,1];
Beq = n;
LB = [2;2;2;2;2;2;2;2];
UB = [5;5;5;5;5;5;5;5];
function [c,ceq] = mycon(x)
    c = [];
    ceq = sin(pi*x); % Compute nonlinear equalities at x.
end
Row = fmincon(@(x) x(1)*x(2)*x(3)*x(4)*x(5)*x(6)*x(7)*x(8)-n,X0,[],[],Aeq,Beq,LB,UB,@mycon);

end