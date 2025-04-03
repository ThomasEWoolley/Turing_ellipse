function csd=G_dmathieuq(kf,m,q,xr,csd,varargin);
%     ===============================================================
%     Purpose: Compute the derivative of  Angular Mathieu functions cem'(q,xr)and sem'(q,xr)  
%  
%     Input :  KF  --- Function code
%     KF=1 for computing even Angular Mathieu function: cem'(q,xr)
%     KF=2 for computing odd Angular Mathieu function: sem'(q,xr)
%     m   --- Order of Mathieu functions
%     q   --- Vector Parameter of Mathieu functions(q can be Positive, Negative or Complex number)
%     x   --- Argument of Mathieu functions(in radian)
%     Output:  csd--- cem'(q,xr)or sem'(q,xr)
%    
%    Editor: Da Ma, Southeast University, China
%     ===============================================================
kd=[];a=[];fg=[];
 lq=length(q);
 fg=zeros(lq,251);
eps=1.0d-14;
if(kf == 1&m == 2.*fix(m./2))kd=1; end;
if(kf == 1&m ~= 2.*fix(m./2))kd=2; end;
if(kf == 2&m ~= 2.*fix(m./2))kd=3; end;
if(kf == 2&m == 2.*fix(m./2))kd=4; end;

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
%====================================================
%computate the devitation of angle mathieu function
%====================================================
csd=zeros(1,lq);
for vm=1:lq
     if imag(q(vm))==0&real(q(vm))<0
        csd(vm)=nan;
     else
       for  k=1:km(vm);
          if(kd == 1);
               csd(vm)=csd(vm)-(2.*k-2).*fg(vm,k).*sin((2.*k-2).*xr);
               elseif(kd == 2);
                    csd(vm)=csd(vm)-(2.*k-1).*fg(vm,k).*sin((2.*k-1).*xr);
                    elseif(kd == 3);
                          csd(vm)=csd(vm)+(2.*k-1).*fg(vm,k).*cos((2.*k-1).*xr);
               elseif(kd == 4);
                    csd(vm)=csd(vm)+2.0d0.*k.*fg(vm,k).*cos(2.*k.*xr);
           end;
      if(k >= ic&abs(fg(vm,k))< abs(csd(vm)).*eps)
          break; 
      end;
   end;
     end
end

return;
end