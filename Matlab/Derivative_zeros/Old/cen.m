%   FUNCTION [Y,LAM]= CEN(X,Q,N) - ELLIPTIC COSINE  v1.0
%   Nth EVEN ANGULAR MATHIEU FUNCTION cen(x,q,n)
%   
%   INPUTS:     -x= vector of values to compute function
%               -q= scalar elliptic parameter value
%               -n= scalar index (n >=0, integer) (the index of the 
%                   eigenfunction' eg. like sin(n*x) )
%   OUTPUTS:    -y= vector of function values
%               -lam= associated eigenvalue
%   The required Angular Mathieu Function is approximated by a trigonometric 
%   fourier expansion. The sub-functions PI_PERIODIC and two_PI_PERIODIC return 
%   the fourier co-efficients (the eigenvectors of the matrix M) and 
%   associated eigenvalue, lambda.
%   Note:   The 'zero-th' eigenfunction can be calculated by using n=0
%           ce0= cen(x,q,0). Since the rows begin counting from n=1 in the
%           program, 
%
%   Based on the formulation given in 
%       L. Chaos-Cador, E. Ley-Koo (2001) 'Mathieu Functions Revisited:
%           matrix evaluation and generating functions'- Revista 
%           Mexicana de Fisica vol.48(1) p 67-75
%
%   -John Terrance 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MAIN FUNCTION CALL
%   default no. of fourier coefficients is 25
function [y, lam] = cen(x,q,n)
    y= zeros(size(x));
    k=n+1;              %shift index n to count matrix index
    if mod(n,2)==0
        [v,mu,idx]= pi_periodic(x,q);
        ncoeffs= size(v,2);
        lam= mu(k);
        y= y + 1/sqrt(2)*v(1,k)*cos(idx(1)*x);
            for j= 2:ncoeffs
                y= y + v(j,k)*cos(idx(j)*x);
            end
    else 
        [v,mu,idx]= two_pi_periodic(x,q);
        ncoeffs= size(v,2);
        lam= mu(n);
            for j= 1:ncoeffs
                y= y + v(j,k)*cos(idx(j)*x);
            end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PI- PERIODIC FUNCTIONS
%   outputs:    -v =fourier expansion coefficients
%               -mu=eigenvalues
function [v,mu,idx]= pi_periodic(x,q);
    idx = 0:24;                         %maximum coefficient index here is 24
    off_diag = q*ones(1,length(idx)-1); off_diag(1)= sqrt(2)*off_diag(1);
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
    idx = 0:24;                         %maximum coefficient index here is 24
    off_diag = q*ones(1,length(idx)-1);
    M= diag((2*idx+1).^2) + diag(off_diag, -1) + diag(off_diag, 1);
    M(1,1)= M(1,1) + q;
%   COMPUTE EIGENVECTORS
    [u,d]   = eig(M);
    [mu,num]= sort(diag(d));
    v= u(:,num);