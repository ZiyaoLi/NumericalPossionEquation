% this is the entrance script of the experiments conducted during the
% project. 8 experiments are conducted in this script. detailed information
% can be found in the corresponding functions called in this script. if
% there is any confusion, feel free to connect the author through email
% [leeeezy@126.com]

% question 1: gauss elimination method demo.
N=64;
eps=[1,1e-3,1e-5];
for i=1:3
    [err,res,res0,t]=GaussElim(N,eps(i));
    ShowResult(1,N,eps(i),t,-1,err,res,res0,-1);
end

% question 2: line g-s iter. + v-cycle + pcg

% hyper-parameters
v1=1;v2=1;
tol=1e-6;
N0=8;

N=[32,64,128,256,512,1024];
eps=[1,1e-1,1e-3,1e-5,1e-7];
for i=1:6
    u0=ones((N(i)-1)*(N(i)-1),1);
    for j=1:5
        [err,res,t,iter]=PCG(N(i),eps(j),u0,tol,N0,0,v1,v2,0);
        ShowResult(2,N(i),eps(j),t,iter,err,res,-1,tol);
    end
end

% question 3: point g-s iter. + v-cycle + pcg

% the same parameters as question 2.

N=256;
eps=[1,1e-3,1e-5];
u0=ones((N-1)*(N-1),1);
for j=1:3
    [err,res,t,iter]=PCG(N,eps(j),u0,tol,N0,1,v1,v2,0);
    ShowResult(3,N,eps(j),t,iter,err,res,-1,tol);
end

% question 4: how good is ordinary cg?

N=[32,64,128,256,512,1024];
eps=[1,1e-1,1e-3,1e-5,1e-7];
for i=1:6
    u0=ones((N(i)-1)*(N(i)-1),1);
    for j=1:5
        [err,res,t,iter]=PCG(N(i),eps(j),u0,tol,N0,3,v1,v2,0);
        ShowResult(5,N(i),eps(j),t,iter,err,res,-1,tol);
    end
end

% question 5: how good is some naive pcgs? such as M=diag(A)?

N=[32,64,128,256,512,1024];
eps=[1,1e-1,1e-3,1e-5,1e-7];
for i=1:6
    u0=ones((N(i)-1)*(N(i)-1),1);
    for j=1:5
        [err,res,t,iter]=PCG(N(i),eps(j),u0,tol,N0,4,v1,v2,0);
        ShowResult(6,N(i),eps(j),t,iter,err,res,-1,tol);
    end
end

% question 6: how about geometry As instead of IAI's?

N=[32,64,128,256,512,1024];
eps=[1,1e-1,1e-3,1e-5,1e-7];
for i=1:6
    u0=ones((N(i)-1)*(N(i)-1),1);
    for j=1:5
        [err,res,t,iter]=PCG(N(i),eps(j),u0,tol,N0,5,v1,v2,0);
        ShowResult(7,N(i),eps(j),t,iter,err,res,-1,tol);
    end
end

% question 7: how good is Matlab's famous '\'?

N=[32,64,128,256,512,1024];
eps=[1,1e-1,1e-3,1e-5,1e-7];
for i=1:6
    u0=ones((N(i)-1)*(N(i)-1),1);
    for j=1:5
        [err,res,t,iter]=PCG(N(i),eps(j),u0,tol,N0,2,v1,v2,0);
        ShowResult(4,N(i),eps(j),t,iter,err,res,-1,tol);
    end
end

% question 8: a plot of line G-S and point G-S

N=256;eps=1e-5;
u0=ones((N-1)*(N-1),1);
k_max=500;
err=zeros(2,k_max);
res=zeros(2,k_max);
[err(1,:),res(1,:),t1,iter1]=PCG(N,eps,u0,tol,N0,0,v1,v2,0);
[err(2,:),res(2,:),t2,iter2]=PCG(N,eps,u0,tol,N0,1,v1,v2,0);
fprintf('Line GS:  t = %3.3f, iter = %3d\n',t1,iter1);
fprintf('Point GS: t = %3.3f, iter = %3d\n',t2,iter2);
plot(err);
plot(res);
