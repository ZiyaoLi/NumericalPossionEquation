function [ U ] = FormingU( N )

% Forming the column vector U (label)
% U(x, y) = sin(pi * x) * sin(pi * y)
% input:
%   N: the number of grids each row.

n=N-1; % to save time repeatedly computing N-1:

U=zeros(n,n);
for i=1:n
    for j=1:n
        U(i,j)=sin(pi*i/N)*sin(pi*j/N);
    end
end
U=U(:);

end
