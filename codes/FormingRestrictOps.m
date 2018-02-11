function [ I ] = FormingRestrictOps( N,N0 )

% this function produces lift and restrict operators (ops). all restrict 
% ops, i.e. I_h^2h, or more descriptively, flat matrices, are produced.
% the corresponding lift ops is 4I'.

% the algorithm looks very stupid, because it calculates the indices and 
% the value vectors directly to form a sparse matrix, since filling the 
% matrix element by element is very time-costing. the algorithm shall
% better be treated as a black-box with correct output, since it is very 
% confusing.

% the algorithm is considered the fastest based on sparse matrices.

% input:
%   N: the initial size (or number of grids per row) of the matrix.
%   N0: the size of the smallest matrix. since a direct method van be
%       applied to solve the inverse of the smallest matrix, there is no
%       need to calculate the restrict / lift ops.
% output:
%   I: the cell of the ops.

I=cell(0);
cnt=0;

while N>N0
    cnt=cnt+1;
    N=N/2; n=N-1; n2=2*N-1;
    % some operations on the indices and values to form a unit, 
    % a unit that replicates all over the full operator.
    unit_i=linspace(1,n,n);
    unit_j=linspace(1,2*n-1,n);
    template_i=repmat(unit_i,3,1);
    template_i=template_i(:)';
    template_j=[unit_j; unit_j+1; unit_j+2];
    template_j=template_j(:)';
    template_v=[ones(1,n);2.*ones(1,n);ones(1,n)];
    template_v=template_v(:)';
    
    % the final indices and value holders:
    i_final=zeros(1,9*n*n);
    j_final=i_final; v_final=i_final;
    
    % filling up the holder:
    for p=1:3
        a=1;
        if p==2
            a=2;
        end
        w=3*n*n;
        
        % filling the (1,2,1) or (2,4,2) sub-matrix indices
        % sub-matrix indices holders:
        i1=zeros(1,w); j1=i1; v1=i1;
        % indices templates:
        i=template_i; j=template_j;
        for r=1:n
            i1((r-1)*(3*n)+1:r*(3*n))=i;
            j1((r-1)*(3*n)+1:r*(3*n))=j;
            v1((r-1)*(3*n)+1:r*(3*n))=template_v*a;
            i=i+n;    % rows indices move n downward
            j=j+2*n2; % columns indices move 2*n2 rightward
        end
        template_j=template_j+n2;  % columns indices move n2 rightward
        i_final((p-1)*w+1:p*w)=i1;
        j_final((p-1)*w+1:p*w)=j1;
        v_final((p-1)*w+1:p*w)=v1;
    end
    
    % finally! forming the matrix!
    I{cnt}=sparse(i_final,j_final,v_final, ... indices and values
                      n*n, n2*n2 ... size of the matrix
                      )./16;
end


end

