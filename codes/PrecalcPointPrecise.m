function [ As,vs ] = PrecalcPointPrecise( A,I,N )
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
%   vs: the critical elements of A at different levels of the v-cycle 
%       process, stored as [a,b,c,d]. detailed information can be found in 
%       function PointGS().

T=length(I);

As=cell(0);
vs=cell(0);

As{1}=A;
vs{1}=[A(1,1),A(2,1),A(N,1),A(N+1,1)];

for i=1:T
    N=N/2;
    As{i+1}=4.*I{i}*As{i}*I{i}';
    vs{i+1}=[As{i+1}(1,1),As{i+1}(2,1),As{i+1}(N,1),As{i+1}(N+1,1)];
end

end

