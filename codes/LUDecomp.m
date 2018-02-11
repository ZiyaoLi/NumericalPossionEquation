function [ LU ] = LUDecomp( A )
% this function calculates an LU decomp. of the given 3-diagonal matrix A,
% returns the coeffs. calculated.

% input:
%   A: the matrix to be decomposed. note that A must be 3-diagonal, or only
%      the diagonal elements are used to calculate.

% output:
%   LU; a cell of 3 elements (l,u1,u2), are the sub-diagonal of calculated
%       L matrix, diagonal and sub-diagonal of the calculated U matrix.
%       such method to store elements is base on the band-width structure
%       of A, and calculated L and U. note that L is a unit triangular
%       matrix and there is no need to store its diagonal elements.

N=length(A);
l=full(diag(A,-1));
u1=full(diag(A));
u2=full(diag(A,1));

for k=1:N-1
    l(k)=l(k)/u1(k);
    u1(k+1)=u1(k+1)-l(k)*u2(k);
end

LU=cell(0);
LU{1}=l;
LU{2}=u1;
LU{3}=u2;

end

