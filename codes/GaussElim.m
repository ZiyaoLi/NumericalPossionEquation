function [ err,res,res0,t ] = GaussElim( N,eps )

% this function uses traditional gauss elimination method to solve the
% equation Au=F given the number of grids each row (N) and epsilon (eps).

% input:
%   N: the number of grids each row.
%   eps: the given eps of the original differential equation:
%        - u_xx - eps * uyy = f(x,y).

% output:
%   err: the l^2 norm of the error: ||U-U^e||_{l^2}.
%   res: the l^2 norm of the residual: ||F-A*U||_{l^2}.
%   res0: the l^2 norm of the residuals from the label: ||F-A*U^e||_{l^2}.
%   t: the time-cost of the calculation.

% calculating A, F an the label of U:
A0=FormingA(N,eps);
f=FormingF(N,eps);
u_true=FormingU(N);

% residual norm generated from r = A * u_true - f
% i.e. the error introduced with discrete differetial op.
% and the unavoidable error of the gauss elimination method.
res0=norm(A0*u_true-f)/N;

A=A0; % for the L-U decomp. method to alter
u=f; % for the triangular equation solution algorithm to alter

tic; % start timing

n=N-1; % to save time repeatedly computing N-1:
L=n*n; % length of the big matrix A, vector u and f.

% L-U decomp. of A; optimized according to the shape of A.
% % First (N-2)*(N-1) rows with at most N-1 elements to eliminate
for k=1:L-n
    % % % low: the lower-bound of the sub-matrix to be adjusted this time
    % % % (or the max i of elements that will change this time)
    low=k+n;
    % % % right: the right-bound of the sub-matrix to be adjusted this time
    % % % (or the max j of elements that will change this time)
    right=min(k+2*n,L);
    A(k+1:low,k)=A(k+1:low,k)./A(k,k);
    A(k+1:low,k+1:right)=A(k+1:low,k+1:right)-A(k+1:low,k)*A(k,k+1:right);
end
% % Last (N-1) rows with at most N-1-k elements to eliminate: no "low" or 
% % "right" optimization can be implemented.
for k=L-n+1:L-1
    A(k+1:L,k)=A(k+1:L,k)/A(k,k);
    A(k+1:L,k+1:L)=A(k+1:L,k+1:L)-A(k+1:L,k)*A(k,k+1:L);
end

% solve Ly=b equation with band-width optimized Algorithm 1.1.1
for j=1:L-1
    % u(j)=u(j)/L(j,j) is unnecessary since L(j,j)=1
    bound=min(j+N,L); % optimized with the band-width information of L
    u(j+1:bound)=u(j+1:bound)-u(j)*A(j+1:bound,j);
end
% solve Ux=y equation with band-width optimized Algorithm 1.1.2
for j=L:-1:2
    u(j)=u(j)/A(j,j);
    bound=max(1,j-N); % optimized with the band-width information of U
    u(bound:j-1)=u(bound:j-1)-u(j)*A(bound:j-1,j);
end
u(1)=u(1)/A(1,1);

t=toc; % stop timing

% % calculate error and residual norms:
err=norm(u-u_true)/N;
res=norm(A0*u-f)/N;

end

