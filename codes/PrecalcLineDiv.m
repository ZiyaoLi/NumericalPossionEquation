function [ As,Ls,Ds ] = PrecalcLineDiv( A,I,N )
% this function does the pre-calculation for a symm. line g-s iter. As of
% different layers of v-cycle, sub-diagonal sub-matrix of As, and the
% diagonal sub-matrix of As are calculated.

% input:
%   A: the original matrix of the problem Au=f.
%   I: the calculated lift / restrict ops from FormingRestricOps().
%   N: the number of grids each row of the original grid.

% output:
%   As: the calculated cells of A at different levels of the v-cycle 
%       process.
%   Ls: the calculated sub-diagonal sub-matrix of A at different levels of
%       the v-cycle process.
%   Ds: the diagonal sub-matrix of A at different levels of the v-cycle 
%        process.

T=length(I);

As=cell(0);
Ls=cell(0);
Ds=cell(0);

As{1}=A;
Ls{1}=As{1}(N:2*N-2,1:N-1);
Ds{1}=As{1}(1:N-1,1:N-1);

for i=1:T
    N=N/2;
    As{i+1}=4.*I{i}*As{i}*I{i}';
    Ls{i+1}=As{i+1}(N:2*N-2,1:N-1);
    Ds{i+1}=As{i+1}(1:N-1,1:N-1);
end

end

