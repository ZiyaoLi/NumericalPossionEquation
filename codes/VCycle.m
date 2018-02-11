function [ u ] = VCycle( As,r,N,method,Is,inv_Ac,const1,const2,v1,v2 )

% this function calculates the eqation Au=r using v-cycle method,
% implementing a line g-s or a point g-s smoother based on 'method'
% argument. some constants are calculated forehead to save time.

% input:
%   As: a cell of the As at different levels.
%   r: the target vector of the equation.
%   N: the size of the smoothest grids.
%   method: which smoother is to choose, 0 for line g-s based on LU decomp,
%           1 for point g-s and 2 for line g-s method based on a '\' op.
%   Is: the pre-calculated restrict ops.
%   inv_Ac: the precalculated direct inverse of the coarsest grids A.
%   const1: the sub-diagonal sub-matrices of different As for line g-s, and
%           critical values for point g-s.
%   const2: the LU decomps. for line g-s with LU, Ds for line g-s with '\'
%           and none for point g-s.
%   v1: the time of iterations forward.
%   v2: the time of iterations backward.

% output:
%   u: the calculated u.

Fs=cell(0);
us=cell(0);

T=length(Is);
Fs{1}=r;
switch method
    case 0
        us{1}=LineGS_LU(zeros((N-1)*(N-1),1),r,N,const1{1},const2{1},v1);
    case 1
        us{1}=PointGS(zeros((N-1)*(N-1),1),r,N,const1{1},v1);
    case 2
        us{1}=LineGS_Div(zeros((N-1)*(N-1),1),r,N,const1{1},const2{1},v1);
    case 5
        us{1}=LineGS_Dir(zeros((N-1)*(N-1),1),r,N,const1{1},v1);
end
r_tmp=Fs{1}-As{1}*us{1};

for i=2:T
    N=N/2;
    Fs{i}=Is{i-1}*r_tmp;
    switch method
        case 0
            us{i}=LineGS_LU(zeros((N-1)*(N-1),1),Fs{i},N,const1{i},const2{i},v1);
        case 1
            us{i}=PointGS(zeros((N-1)*(N-1),1),Fs{i},N,const1{i},v1);
        case 2
            us{i}=LineGS_Div(zeros((N-1)*(N-1),1),Fs{i},N,const1{i},const2{i},v1);
        case 5
            us{i}=LineGS_Dir(zeros((N-1)*(N-1),1),Fs{i},N,const1{i},v1);
    end
    r_tmp=Fs{i}-As{i}*us{i};
end
Fs{T+1}=Is{T}*r_tmp;
us{T+1}=inv_Ac*Fs{T+1};
for i=T:-1:1
    u=us{i}+(4.*Is{i}')*us{i+1};
    switch method
        case 0
            us{i}=LineGS_LU(u,Fs{i},N,const1{i},const2{i},v2);
        case 1
            us{i}=PointGS(u,Fs{i},N,const1{i},v2);
        case 2
            us{i}=LineGS_Div(u,Fs{i},N,const1{i},const2{i},v2);
        case 5
            us{i}=LineGS_Dir(u,Fs{i},N,const1{i},v2);
    end
    N=N*2;
end

% retrieve the final output
u=us{1};

end

