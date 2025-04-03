function d4r=dmmathieu4q(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute the derivative of radial Mathieu functions of the forrth kind Mcm'(q,x)and Msm'(q,x)
%
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm'(q,x)
%     KF=2 for computing odd function Msm'(q,x)
%     m  --- Order of Mathieu functions
%     q  ---Vector Parameter of Mathieu functions(q > 0)
%     x  ---Arguments of rsdial Mathieu functions
%     Output:  f4r --- Mcm'(q,x)or Msm'(q,x)
%
% Editor: Da Ma, Southeast University, China
%     ==============================================================
d1r=dmmathieu1q(kf,m,q,x,varargin);
d2r=dmmathieu2q(kf,m,q,x,varargin);
d4r=d1r-i.*d2r;
end