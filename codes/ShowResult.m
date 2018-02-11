function [  ] = ShowResult( type,N,eps,time,iter,err,res,res0,tol )
% this functions show the calculation results according to the inputs.
% the meaning of inputs are clearly shown through variable names.

switch type
    case 1  % gauss elimination
        fprintf('\n');
        fprintf('### Gauss Elimination Method ###\n');
        fprintf('# N = %d, eps = %3.0e\n',N,eps);
        fprintf('# cpu_time = %3.4fs\n',time);
        fprintf('# error_norm = %3.3e, residual_norm = %3.3e\n',err,res);
        fprintf('# residual_norm0 = %3.3e\n',res0);
        fprintf('################################\n\n');
    case 2  % line g-s (LU) + v-cycle + pcg
        fprintf('\n');
        fprintf('### PCG with V-Cycle and Line G-S (LU) ###\n');
        fprintf('# N = %d, eps = %3.0e, tol = %3.0e\n',N,eps,tol);
        fprintf('# cpu_time = %3.4fs, iter = %3d\n',time,iter);
        fprintf('##########################################\n\n');
    case 3 % point g-s + v-cycle + pcg
        fprintf('\n');
        fprintf('### PCG with V-Cycle and Point G-S ###\n');
        fprintf('# N = %d, eps = %3.0e, tol = %3.0e\n',N,eps,tol);
        fprintf('# cpu_time = %3.4fs, iter = %3d\n',time,iter);
        fprintf('######################################\n\n');
    case 4 % line g-s (Div) + v-cycle + pcg
        fprintf('\n');
        fprintf('### PCG with V-Cycle and Line G-S (Div) ###\n');
        fprintf('# N = %d, eps = %3.0e, tol = %3.0e\n',N,eps,tol);
        fprintf('# cpu_time = %3.4fs, iter = %3d\n',time,iter);
        fprintf('###########################################\n\n');
    case 5 % ordinary CG
        fprintf('\n');
        fprintf('### Ordinary CG Method ###\n');
        fprintf('# N = %d, eps = %3.0e, tol = %3.0e\n',N,eps,tol);
        fprintf('# cpu_time = %3.4fs, iter = %3d\n',time,iter);
        fprintf('##########################\n\n');
    case 6 % naive PCG
        fprintf('\n');
        fprintf('### Naive PCG Method: Diagonal Matrix Approx. ###\n');
        fprintf('# N = %d, eps = %3.0e, tol = %3.0e\n',N,eps,tol);
        fprintf('# cpu_time = %3.4fs, iter = %3d\n',time,iter);
        fprintf('#################################################\n\n');
    case 7 % line g-s (Div) + v-cycle + pcg
        fprintf('\n');
        fprintf('### PCG with V-Cycle and Line G-S (Dir) ###\n');
        fprintf('# N = %d, eps = %3.0e, tol = %3.0e\n',N,eps,tol);
        fprintf('# cpu_time = %3.4fs, iter = %3d\n',time,iter);
        fprintf('###########################################\n\n');
end


end

