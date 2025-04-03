function f2r=G_mmathieu2q(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute modified Mathieu functions of the second kinds Mcm(q,x)and Msm(q,x)
%     
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm(x,q)
%     KF=2 for computing odd function Msm(x,q)
%     m  --- Order of Mathieu functions
%     q  --- Vector Parameter of Mathieu functions(q can be Positive, Negative or Complex number)
%     x  ---Argument of Mathieu functions
%     Output:   F2R --- Mcm(x,q)or Msm(x,q)
%
%    Editor: Da Ma, Southeast University, China.
%     ==============================================================
kd=[];a=[];fg=[];km=[];u1=[];nm=[];bj1=[];dj1=[];by1=[];dy1=[];u2=[];bj2=[];dj2=[];by2=[];dy2=[];
nr=length(x);
lq=length(q);
 fg=zeros(lq,251);
 a=0;
bj1=zeros(1,lq);
 bj2=zeros(1,lq);
 by1=zeros(1,251+1);
 by2=zeros(1,251+1);
c1=0;
 c2=0;
 dj1=zeros(1,251+1);
 dj2=zeros(1,251+1);
 dy1=zeros(1,251+1);
dy2=zeros(1,251+1);
 eps=0;
 qm=0;
 u1=0;
 u2=0;
 w1=0;
 w2=0;
 ic=0;
 k=0;
 kd=0;
 km=0;
 nm=0;
eps=1.0d-14;
if(kf==1 & m==2.*fix(m./2));
kd=1;
end;
if(kf==1 & m~=2.*fix(m./2));
kd=2;
end;
if(kf==2 & m~=2.*fix(m./2));
kd=3;
end;
if(kf==2 & m==2.*fix(m./2));
kd=4;
end;
for k=1:lq
 if(abs(q(k)) <= 1.0d0);
     qm(k)=7.5+56.1.*sqrt(abs(q(k)))-134.7.*abs(q(k))+90.7.*sqrt(abs(q(k))).*abs(q(k));
       else;
         qm(k)=17.0+3.1.*sqrt(abs(q(k)))-.126.*abs(q(k))+.0037.*sqrt(abs(q(k))).*abs(q(k));
   end;
   km(k)=fix(qm(k)+0.5.*fix(m));
end
for k=1:lq
fc=fcoef1(kd,fix(m),q(k)).';
fg([k],1:length(fc))=fc;
end
ic=fix(fix(fix(m)./2)+1);
if(kd==4);
ic=fix(fix(m)./2);
end;
kr=1:lq;
c1=exp(-x);
c2=exp(x);
u1(kr)=sqrt(q(kr)).*c1;
u2(kr)=sqrt(q(kr)).*c2;
for kr=1:lq 
k=1:km+2;
bj1(k,kr)=besselj(fix(k)-1,u1(kr));
by2(k,kr)=bessely(fix(k)-1,u2(kr));
end;
%======================================================================
%calculation the second kind modified mathieu functions
%======================================================================
f2r=zeros(1,lq);
for kr=1:lq;
for k=1 : km;
if(kd==1);
f2r(kr)=f2r(kr)+(-1).^(fix(ic)+fix(k)).*fg(kr,fix(k)).*bj1(fix(k)-1+1,kr).*by2(fix(k)-1+1,kr);
elseif(kd==2 | kd==3);
f2r(kr)=f2r(kr)+(-1).^(fix(ic)+fix(k)).*fg(kr,fix(k)).*(bj1(fix(k)-1+1,kr).*by2(fix(k)+1,kr)+(-1).^fix(kd).*bj1(fix(k)+1,kr).*by2(fix(k)-1+1,kr));
else;
f2r(kr)=f2r(kr)+(-1).^(fix(ic)+fix(k)).*fg(kr,fix(k)).*(bj1(fix(k)-1+1,kr).*by2(fix(k)+1+1,kr)-bj1(fix(k)+1+1,kr).*by2(fix(k)-1+1,kr));
end;
if(k>=5 & abs(f2r(kr)-w1)<abs(f2r(kr)).*eps);
break;
end;
w1=f2r(kr);
end;
f2r(kr)=f2r(kr)./fg(kr,1);
end;
end  