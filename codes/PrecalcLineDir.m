function [ As,LUs ] = PrecalcLineDir( I,N,eps )
% this function does the pre-calculation for a symm. line g-s iter. As of
% different layers of v-cycle, sub-diagonal sub-matrix of As, and the
% diagonal sub-matrix of As are calculated.

% Note that As are not formed according to the restrict / lift ops, but a
% corresponding differential matrix to N.

% input:
%   I: the calculated lift / restrict ops from FormingRestricOps().
%   N: the number of grids each row of the original grid.
%   eps: the eps shown in the equation.

% output:
%   As: the calculated cells of A at different levels of the v-cycle 
%       process.
%   LUs: the LU decomps. of the diagonal sub-matrix of A at different
%        levels of the v-cycle process.

T=length(I);

As=cell(0);
LUs=cell(0);

for i=1:T+1
    As{i}=FormingA(N,eps);
    LUs{i}=LUDecomp(As{i}(1:(N-1),1:(N-1)));
    N=N/2;
end

end

