function [ As,Ls,LUs ] = PrecalcLineLU( A,I,N )
% this function does the pre-calculation for a symm. line g-s iter. As of
% different layers of v-cycle, sub-diagonal sub-matrix of As, and the LU
% decomps. of the diagonal sub-matrix of As are calculated.

% input:
%   A: the original matrix of the problem Au=f.
%   I: the calculated lift / restrict ops from FormingRestricOps().
%   N: the number of grids each row of the original grid.

% output:
%   As: the calculated cells of A at different levels of the v-cycle 
%       process.
%   Ls: the calculated sub-diagonal sub-matrix of A at different levels of
%       the v-cycle process.
%   LUs: the LU decomps. of the diagonal sub-matrix of A at different
%        levels of the v-cycle process.

T=length(I);

As=cell(0);
Ls=cell(0);
LUs=cell(0);

As{1}=A;
LUs{1}=LUDecomp(A(1:N-1,1:N-1));
Ls{1}=As{1}(N:2*N-2,1:N-1);

for i=1:T
    N=N/2; n=N-1;
    As{i+1}=4.*I{i}*As{i}*I{i}';
    LUs{i+1}=LUDecomp(As{i+1}(1:n,1:n));
    Ls{i+1}=As{i+1}(n+1:2*n,1:n);
end

end

