%   FUNCTION [Y,LAM]= CEN(X,Q,N) - ELLIPTIC SINE  v1.0
%   Nth ODD ANGULAR MATHIEU FUNCTION sen(x,q,n)
%   
%   INPUTS:     -x= vector of values to compute function
%               -q= scalar elliptic parameter value
%               -n= positive (n>0) scalar index (integer) (the index of the 
%                   eigenfunction' eg. like sin(n*x) )
%   OUTPUTS:    -y= vector of function values
%               -lam= associated eigenvalue
%   The required Angular Mathieu Function is approximated by a trigonometric 
%   fourier expansion. The sub-functions PI_PERIODIC and two_PI_PERIODIC return 
%   the fourier co-efficients (the eigenvectors of the matrix M) and 
%   associated eigenvalue, lambda.
%
%   Note:   There is no 'zero-th' odd mathieu function as there is in the
%           even (ce) case.
%
%   Based on the formulation given in 
%       L. Chaos-Cador, E. Ley-Koo (2001) 'Mathieu Functions Revisited:
%           matrix evaluation and generating functions'- Revista 
%           Mexicana de Fisica vol.48(1) p 67-75
%
%   -John Terrance 2007mathieu_example.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MAIN FUNCTION CALL
%   default no. of fourier coefficients is 25
function [y, lam] = sen(x,q,n)
y= zeros(size(x));
if mod(n,2)==0
    [v,mu,idx]= pi_periodic(x,q);
    ncoeffs= size(v,2);
    lam= mu(n);
        for j= 1:ncoeffs
            y= y + v(j,n)*sin(idx(j)*x);
        end
 else
    [v,mu,idx]= two_pi_periodic(x,q);
    ncoeffs= size(v,2);
    lam= mu(n);
        for j= 1:ncoeffs
            y= y + v(j,n)*sin(idx(j)*x);
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PI-PERIODIC FUNCTIONS
%   outputs:    -v =fourier expansion coefficients
%               -mu=eigenvalues
function [v,mu,idx]= pi_periodic(x,q);
idx = 1:25;                             %maximum coefficient index here is 24
off_diag = q*ones(1,length(idx)-1);     
M= diag((2*idx).^2) + diag(off_diag, -1) + diag(off_diag, 1);
%   COMPUTE EIGENVECTORS
[u,d]   = eig(M);
[mu,num]= sort(diag(d));
v= u(:,num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2*PI- PERIODIC FUNCTIONS
%   outputs:    -v =fourier expansion coefficients
%               -mu=eigenvalues
function [v,mu,idx]= two_pi_periodic(x,q);
idx = 0:24;                             %maximum coefficient index here is 24
off_diag = q*ones(1,length(idx)-1);
M= diag((2*idx+1).^2) + diag(off_diag, -1) + diag(off_diag, 1);
M(1,1)= M(1,1) - q;
%   COMPUTE EIGENVECTORS
[u,d]   = eig(M);
[mu,num]= sort(diag(d));
v= u(:,num);