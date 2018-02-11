function [ b ] = LUSolve( LU,b )

% this function computes the solution of the equation LUx=b, where
%   |  1  0  0  0  0 |
%   | l1  1  0  0  0 |
% L=|  0 l2  1  0  0 |, 
%   |  0  0 l3  1  0 |
%   |  0  0  0 l4  1 |
% 
%   | u11 u21  0   0   0  |
%   |  0  u12 u22  0   0  |
% U=|  0   0  u13 u23  0  |
%   |  0   0   0  u14 u24 |
%   |  0   0   0   0  u15 |
% and the inputs l, u1, u2 are elements in L and U, namely:
%        l=(l1,l2,l3,l4)', 
%       u1=(u11,u12,u13,u14,u15)', 
% and   u2=(u21,u22,u23,u24)'.

% input:
%   LU: the calculated LU decomp. from LUDecomp(). for more detail, please
%       check the corresponding function.
%   b: the target vector of the equation.

% output:
%   b: the calculated vector, stored in b.

% retrieve coefficients from the cell
l=LU{1};
u1=LU{2};
u2=LU{3};
n=length(u1);

% solve Ly=b equation with 3-diagonal optimized Algorithm 1.1.1
for j=1:n-1
    b(j+1)=b(j+1)-b(j)*l(j);
end

% solve Ux=y equation with band-width optimized Algorithm 1.1.2
for j=n:-1:2
    b(j)=b(j)/u1(j);
    b(j-1)=b(j-1)-b(j)*u2(j-1);
end
b(1)=b(1)/u1(1);

end
