function A = outerproduct(a,b,c,d,e)
i = length(a);
j= length(b);
k = length(c);
% check input requirements
n = nargin;
switch n
    case 2
        A = kron(a,b');
    case 3
        A = reshape(kron(c',kron(a,b')),[i,j,k]);
    case 4        
        l = length(d);
        A = reshape(kron(d',kron(c',kron(a,b'))),[i,j,k,l]);
    case 5
        l = length(d);
        m = length(e);
        A = reshape(kron(e',kron(d',kron(c',kron(a,b')))),[i,j,k,l,m]);
end

end