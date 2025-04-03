function f3r=G_mmathieu3q(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute modified Mathieu functions of the third kinds Mcm(q,x)and Msm(q,x)
%     
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm(x,q)
%     KF=2 for computing odd function Msm(x,q)
%     m  --- Order of Mathieu functions
%     q  --- Vector Parameter of Mathieu functions(q can be Positive, Negative or Complex number)
%     x  ---Argument of Mathieu functions
%     Output:   F3R --- Mcm(x,q)or Msm(x,q)
%
%    Editor: Da Ma, Southeast University, China.
%     ==============================================================

f1r=G_mmathieu1q(kf,m,q,x,varargin);
f2r=G_mmathieu2q(kf,m,q,x,varargin);
f3r=f1r+i.*f2r;
end