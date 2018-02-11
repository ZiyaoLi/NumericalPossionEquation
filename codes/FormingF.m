function [ F ] = FormingF( N, eps )

% Forming the column vector F
% F(x, y) = (1 + eps) * pi ^ 2 * sin(pi * x) * sin(pi * y)

n=N-1; % to save time repeatedly computing N-1:

F=zeros(n,n);
for i=1:n
    for j=1:n
        F(i,j)=(1+eps)*pi*pi*sin(pi*i/N)*sin(pi*j/N);
    end
end
F=F(:);

end