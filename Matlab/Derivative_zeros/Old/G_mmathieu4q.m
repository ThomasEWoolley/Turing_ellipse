function f4r=G_mmathieu4q(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute radial Mathieu functions of the fourth kind Mcm(q,x)and Msm(q,x)
%
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm'(q,x)
%     KF=2 for computing odd function Msm'(q,x)
%     m  --- Order of Mathieu functions
%     q  ---Vector Parameter of Mathieu functions(q can be Positive, Negative or Complex number)
%     x  ---Arguments of rsdial Mathieu functions
%     Output:  f4r --- Mcm(q,x)or Msm(q,x)
%
% Editor: Da Ma, Southeast University, China
%     ==============================================================
f1r=G_mmathieu1q(kf,m,q,x,varargin);
f2r=G_mmathieu2q(kf,m,q,x,varargin);
f4r=f1r-i.*f2r;
end