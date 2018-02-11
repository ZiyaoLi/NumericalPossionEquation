function [ u ] = LineGS_Div( u, r, N, L, D, iter )
% using iter-step symmetric line gauss-seidel iterative method
% as the smoother of v-cycle precoditioner solving "Au=r",
% where
%           | D L 0 0 0 |
%           | L D L 0 0 |
%         A=| 0 L D L 0 |
%           | 0 0 L D L |
%           | 0 0 0 L D |
% and D, L are (N-1)*(N-1) square matrices. using matlab's fastest solver 
% of linear equations: '\'. no precalculation (pre LU decomp.) implemented.

% input:
%   u: the initial value of u.
%   r: the target vector of Au=r.
%   N: the current number of grids per row, or the number of rows of D +1.
%   L: the sub-diagonal sub-matrix of A.
%   D: the diagonal sub-matrix of A.
%   iter: the number of iteration.

% output
%   u: the calculated u.

% to save time of massive calculations of N-1) in the cycle:
n=N-1;

for t=1:iter
    
    % order iteration:
    % special case in the first row
    u(1:n)=D\(r(1:n)-L*u(n+1:2*n));
    
    for i=2:n-1
        s=(i-1)*n;
        u(s+1:s+n)=D\(r(s+1:s+n)-L*(u(s-n+1:s)+u(s+n+1:s+2*n)));
    end
    
    % special case in the last row
    s=(n-1)*n;
    u(s+1:s+n)=D\(r(s+1:s+n)-L*u(s-n+1:s));
    
    % reverse-order iteration:
    % very similar to the order iteration.
    u(s+1:s+n)=D\(r(s+1:s+n)-L*u(s-n+1:s));
    for i=n-1:-1:2
        s=(i-1)*n;
        u(s+1:s+n)=D\(r(s+1:s+n)-L*(u(s-n+1:s)+u(s+n+1:s+2*n)));
    end
    u(1:n)=D\(r(1:n)-L*u(n+1:2*n));
end

end

