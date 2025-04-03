function csf=G_mathieuq(kf,m,q,xr,csf,varargin)
%     ===============================================================
%     Purpose: Compute Angular Mathieu functions cem(q,xr)and sem(q,xr)
%  
%     Input :  KF  --- Function code
%     KF=1 for computing even Angular Mathieu function: cem(q,xr)
%     KF=2 for computing odd Angular Mathieu function: sem(q,xr)
%     m   --- Order of Mathieu functions
%     q   --- Vector Parameter of Mathieu functions(q can be Positive, Negative or Complex number)
%     x   --- Argument of Mathieu functions(in radian)
%     Output:  CSF --- cem(q,xr)or sem(q,xr)
%    
%Editor: Da Ma,Southeast University, China
%     ===============================================================
kd=[];a=[];fg=[];
lq=length(q);
 fg=zeros(lq,251);
eps=1.0d-14;
if(kf == 1&m == 2.*fix(m./2))kd=1; end;
if(kf == 1&m ~= 2.*fix(m./2))kd=2; end;
if(kf == 2&m ~= 2.*fix(m./2))kd=3; end;
if(kf == 2&m == 2.*fix(m./2))kd=4; end;
%%========================================================================
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
ic=fix(fix(m)./2)+1;
csf=zeros(1,lq);
for vm=1:lq
    if imag(q(vm))==0&real(q(vm))<0
        csf(vm)=Mathieu_Q(kf,m,q(vm),xr);
    else
 for  k=1:km(vm);
   if(kd == 1);
       csf(vm)=csf(vm)+fg(vm,k).*cos((2.*k-2).*xr);
   elseif(kd == 2);
       csf(vm)=csf(vm)+fg(vm,k).*cos((2.*k-1).*xr);
    elseif(kd == 3);
       csf(vm)=csf(vm)+fg(vm,k).*sin((2.*k-1).*xr);
     elseif(kd == 4);
       csf(vm)=csf(vm)+fg(vm,k).*sin(2.*k.*xr);
     end;
   if(k >= ic&abs(fg(vm,k))< abs(csf(vm)).*eps)
       break; 
   end;
 end;
    end 
end;
return;
end
