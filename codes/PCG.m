function [ err,res,t,iter ] = PCG( N,eps,u0,tol,N0,method,v1,v2,verbose )

% this function uses preconditioning conjugate gradient (PCG) method to 
% solve the equation Au=F given the number of grids each row (N) and 
% epsilon (eps). v-cycle is used as the preconditioner, and line g-s iter.
% is used as the smoother.

% input:
%   N: the number of grids each row.
%   eps: the given eps of the original differential equation:
%        - u_xx - eps * uyy = f(x,y).
%   u0: the initial value given.
%   tol: tolerant reletive residual or standard of halting, i.e. halt when 
%        ||r|| / ||r0|| < tol.
%   N0: the size of the coarsest v-cycle grids; the u is calculated
%       directly from the inverse of the matrix.
%   method: whether a line g-s or a point g-s is implemented. 0 for line;
%           1 for point; 2 for line using '\' op instead of LU decomp; 
%           3 for ordinary cg; 4 for preconditioning with D=diag(A).
%   v1, v2: the time of iterations forward and backward in v-cycle.
%   verbose: whether print the detailed information of each iteration.

% output:
%   err: the history of l^2 norms of the error: ||U-U^e||_{l^2}.
%   res: the history of l^2 norms of the residual: ||F-A*U||_{l^2}.
%   t: the time-cost of the calculation.
%   iter: number of the iterations of PCG.

% calculating A, F an the label of U:
A=FormingA(N,eps);
f=FormingF(N,eps);
u_true=FormingU(N);

r0=f-A*u0;

tic; % start timing

% get the lift/restrict ops beforehead to save time.
Is=FormingRestrictOps(N,N0);

% pre-calculation of some constants needed in the v-cycle process.
% for detailed information, please check the corresponding pre-calculation
% function descriptions or the v-cycle function descriptions.
switch method
    case 0
        [As,Ls,LUs]=PrecalcLineLU(A,Is,N);
        inv_Ac=full(As{length(As)})^-1;
        k_max=500;
    case 1
        [As,vs]=PrecalcPointPrecise(A,Is,N);
        inv_Ac=full(As{length(As)})^-1;
        k_max=500;
    case 2
        [As,Ls,Ds]=PrecalcLineDiv(A,Is,N);
        inv_Ac=full(As{length(As)})^-1;
        k_max=500;
    case 3
        k_max=10000;
    case 4
        n=N-1;
        D=A(1:n,1:n);
        LU=LUDecomp(D);
        k_max=10000;
    case 5
        [As,LUs]=PrecalcLineDir(Is,N,eps);
        inv_Ac=full(As{length(As)})^-1;
        k_max=500;
end

C=norm(r0)*tol;  % the absolute halting standard.

% vehicles for recording error and residual norms:
err=zeros(k_max,1); res=err;

% PCG process.
% mainly identical with Algorithm 5.4.1, and the part solving Mz=r is 
% replaced with the v-cycle process.
u=u0; k=0; r=r0; % initialization of the pcg method.

while(norm(r)>C && k<k_max)
    if verbose
        % a glimpse on the current iter and the current residual norm rate.
        fprintf('iter %3d: current_norm/halt_norm = %3.2f\n\n',k,norm(r)/C)
    end
    k=k+1; 
    
    % recording the norm of error and residual during the iteration.
    err(k)=norm(u-u_true)/N;
    res(k)=norm(r)/N;
    
    % v-cycle preconditioning.
    % detailed information to be found in the corresponding func.
    switch method
        case 0
            z=VCycle(As,r,N,method,Is,inv_Ac,Ls,LUs,v1,v2);
        case 1
            z=VCycle(As,r,N,method,Is,inv_Ac,vs,0,v1,v2);
        case 2
            z=VCycle(As,r,N,method,Is,inv_Ac,Ls,Ds,v1,v2);
        case 3
            z=r;
        case 4
            z=zeros(n*n,1);
            for i=1:n
                z((i-1)*n+1:i*n)=LUSolve(LU,r((i-1)*n+1:i*n));
            end
        case 5
            z=VCycle(As,r,N,method,Is,inv_Ac,LUs,0,v1,v2);
    end
    
    if k==1
        p=z; rho=r'*z;
    else
        rho_tilde=rho; rho=r'*z;
        beta=rho/rho_tilde; p=z+beta.*p;
    end
    
    w=A*p; alpha=rho/(p'*w);
    u=u+alpha.*p; r=r-alpha.*w;
end

t=toc; % stop timing

iter=k;

end

