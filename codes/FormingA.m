function [ A ] = FormingA( N,eps )

% this function returns the A matrix corresponding to given n and eps.
% input:
%   N: the number of grids each row.
%   eps: the given eps of the original differential equation:
%        - u_xx - eps * uyy = f(x,y).

% to save time repeatedly computing N-1:
n=N-1;

% form a unit diagonal matrix I_n
I_n=speye(n);

% form a unit sub-diagonal matrix I_nsub
i_sub=linspace(1,n-1,n-1);
j_sub=i_sub+1;
I_nsub=sparse(i_sub,j_sub,ones(1,n-1),n,n);

% form diagonal sub-matrices a_n:
a_n=2*(1+eps).*I_n-eps.*I_nsub-eps.*I_nsub';

% form the final A:
A=kron(I_n,a_n)      ...  % diagonal sub-matrices
    -kron(I_nsub,I_n)...  % upper sub-diagonal sub-matrices
    -kron(I_nsub',I_n);   % lower sub-diagonal sub-matrices

A=A.*(N*N);

end

