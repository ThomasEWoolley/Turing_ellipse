function f4r=mmathieu4q(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute radial Mathieu functions of the fourth kind Mcm(q,x)and Msm(q,x)
%
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm'(q,x)
%     KF=2 for computing odd function Msm'(q,x)
%     m  --- Order of Mathieu functions
%     q  ---Vector Parameter of Mathieu functions(q > 0)
%     x  ---Arguments of rsdial Mathieu functions
%     Output:  f4r --- Mcm(q,x)or Msm(q,x)
%
% Editor: Da Ma, Southeast University, China
%     ==============================================================
f1r=mmathieu1q(kf,m,q,x,varargin);
f2r=mmathieu2q(kf,m,q,x,varargin);
f4r=f1r-i.*f2r;
end