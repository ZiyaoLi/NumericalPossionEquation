function [ u ] = PointGS( u, r, N, v, iter )
% using iter-step symmetric point gauss-seidel iterative method
% as the smoother of v-cycle precoditioner solving "Au=r",
% where
%           | D L 0 0 0 |
%           | L D L 0 0 |
%         A=| 0 L D L 0 |               ,
%           | 0 0 L D L |
%           | 0 0 0 L D |(N-1)*(N-1)^2
%
%           | a b 0 0 0 |
%           | b a b 0 0 |
%         D=| 0 b a b 0 |               and
%           | 0 0 b a b |
%           | 0 0 0 b a |(N-1)^2
%
%           | c d 0 0 0 |
%           | d c d 0 0 |
%         L=| 0 d c d 0 |               .
%           | 0 0 d c d |
%           | 0 0 0 d c |(N-1)^2
% output the smoothed 'u'.

% the complexity of this algorithm is O(8*iter*N^2), or linear to the 
% number of rows of A and the elements of u.

% input:
%   u: the initial value of u.
%   r: the target vector of Au=r.
%   N: the current number of grids per row, or the number of rows of D +1.
%   v: the critical elements in D and L: v=[a,b,c,d];
%   iter: the number of iteration.

% output:
%   u: the smoothed vector of Au=r.

% to save time of massive additions of N+(-1) in the cycle:
n=N-1;

a=v(1);b=v(2);c=v(3);d=v(4);

for t=1:iter
    
    % order iteration:
    % manually figures the position of non-zero elements in A to save calculating time.
    u(1)=(r(1)-b*u(2)-c*u(1+n)-d*u(2+n))/a;  % 3 non-zero elements of A in this row.
    for j=2:n-1
        u(j)=(r(j)-b*(u(j-1)+u(j+1))-c*u(j+n)-d*(u(j+n-1)+u(j+n+1)))/a;  % 5 non-zero elements of A in this row.
    end
    u(n)=(r(n)-b*u(n-1)-c*u(n+n)-d*u(n+n-1))/a;  % 3 non-zero elements of A in this row.
    for i=2:n-1
        s=(i-1)*n;
        u(s+1)=(r(s+1)-b*u(s+2)-c*(u(s+1-n)+u(s+1+n))-d*(u(s+1-n+1)+u(s+1+n+1)))/a;  % 5 non-zero elements of A in this row.
        for j=2:n-1  % 8 non-zero elements of A in this row.
            u(s+j)=(r(s+j)-b*(u(s+j+1)+u(s+j-1))-c*(u(s+j-n)+u(s+j+n))-d*(u(s+j-n-1)+u(s+j-n+1)+u(s+j+n-1)+u(s+j+n+1)))/a;
        end
        u(s+n)=(r(s+n)-b*u(s+n-1)-c*(u(s+n-n)+u(s+n+n))-d*(u(s+n-n-1)+u(s+n+n-1)))/a;  % 5 non-zero elements of A in this row.
    end
    s=(n-1)*n;
    u(s+1)=(r(s+1)-b*u(s+1+1)-c*u(s+1-n)-d*u(s+1-n+1))/a;  % 3 non-zero elements of A in this row.
    for j=2:n-1
        u(s+j)=(r(s+j)-b*(u(s+j-1)+u(s+j+1))-c*u(s+j-n)-d*(u(s+j-n-1)+u(s+j-n+1)))/a;  % 5 non-zero elements of A in this row.
    end
    u(s+n)=(r(s+n)-b*u(s+n-1)-c*u(s)-d*u(s-1))/a;  % 3 non-zero elements of A in this row.
    
    % reverse-order iteration:
    % very similar to the order iteration.
    u(s+n)=(r(s+n)-b*u(s+n-1)-c*u(s)-d*u(s-1))/a;  % 3
    for j=n-1:-1:2
        u(s+j)=(r(s+j)-b*(u(s+j-1)+u(s+j+1))-c*u(s+j-n)-d*(u(s+j-n-1)+u(s+j-n+1)))/a;  % 5
    end
    u(s+1)=(r(s+1)-b*u(s+2)-c*u(s-n+1)-d*u(s-n+2))/a;  % 3
    for i=n-1:-1:2
        s=(i-1)*n;
        u(s+n)=(r(s+n)-b*u(s+n-1)-c*(u(s)+u(s+2*n))-d*(u(s-1)+u(s+2*n-1)))/a;  % 5
        for j=n-1:-1:2  % 8
            u(s+j)=(r(s+j)-b*(u(s+j+1)+u(s+j-1))-c*(u(s+j-n)+u(s+j+n))-d*(u(s+j-n-1)+u(s+j-n+1)+u(s+j+n-1)+u(s+j+n+1)))/a;
        end
        u(s+1)=(r(s+1)-b*u(s+2)-c*(u(s-n+1)+u(s+n+1))-d*(u(s-n+2)+u(s+n+2)))/a;  % 5
    end
    u(n)=(r(n)-b*u(n-1)-c*u(2*n)-d*u(2*n-1))/a;  % 3
    for j=n-1:-1:2
        u(j)=(r(j)-b*(u(j-1)+u(j+1))-c*u(j+n)-d*(u(j+n-1)+u(j+n+1)))/a;  % 5
    end
    u(1)=(r(1)-b*u(2)-c*u(n+1)-d*u(n+2))/a;  % 3
end
end

