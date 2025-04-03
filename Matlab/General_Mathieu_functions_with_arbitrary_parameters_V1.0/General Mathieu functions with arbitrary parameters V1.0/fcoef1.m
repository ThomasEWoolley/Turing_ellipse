function fc=fcoef1(kd,m,q,fc,varargin);
%     =====================================================
%     Purpose: Compute expansion coefficients for Angular and Radial Mathieu functions 
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions (q can be Positive, Negative or Complex number)
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,...)
%     KD=2 for cem(x,q,m = 1,3,5,...)
%     KD=3 for sem(x,q,m = 1,3,5,...)
%     KD=4 for sem(x,q,m = 2,4,6,...)
%     A  --- Characteristic value of Mathieu functions for given m and q
%     Output:  FC(k)--- Expansion coefficients of Mathieu
%     functions(k= 1,2,...,KM)
%     FC(1),FC(2),FC(3),... correspond to
%     A0,A2,A4,... for KD=1 case, A1,A3,
%     A5,... for KD=2 case, B1,B3,B5,...
%     for KD=3 case and B2,B4,B6,... for
%     KD=4 case
%
%     Editor: Da Ma, Southeast University, China
%     =====================================================
    if(abs(q) <= 1.0d0);
    qm=7.5+56.1.*sqrt(abs(q))-134.7.*abs(q)+90.7.*sqrt(abs(q)).*abs(q);
    else;
    qm=17.0+3.1.*sqrt(abs(q))-.126.*abs(q)+.0037.*sqrt(abs(q)).*abs(q);
    end;
    km=fix(qm+0.5.*fix(m));
    n=km;
    if kd==1
     e1=q.*ones(n,1);e2=zeros(n,1);e3=q.*ones(n,1);
      for n=1:km
       e2(n,1)=(2*(n-1))^2;
      end
      s3=spdiags([e1,e2,e3],[-1,0,1],n,n);
      s3(1,1)=eps;
      s3(2,1)=2*q;
      s3;
     [v,d]=eig(full(s3));

  fc=v(:,(m+2)/2)./sqrt(1+abs(v(1,(m+2)/2))^2);
elseif kd==2
   e1=q.*ones(n,1);e2=zeros(n,1);e3=q.*ones(n,1);
   for n=1:km
       e2(n,1)=(2*n-1)^2;
   end
   s3=spdiags([e1,e2,e3],[-1,0,1],n,n);
  s3(1,1)=1+q;
  [v,d]=eig(full(s3));
 fc= v(:,(m+1)/2);
elseif kd==3   
   e1=q.*ones(n,1);e2=zeros(n,1);e3=q.*ones(n,1);
   for n=1:km
       e2(n,1)=(2*n-1)^2;
   end
   s3=spdiags([e1,e2,e3],[-1,0,1],n,n);
  s3(1,1)=1-q;
  [v,d]=eig(full(s3));
fc=v(:,(m+1)/2);
elseif kd==4
   e1=q.*ones(n,1);e2=zeros(n,1);e3=q.*ones(n,1);
   for n=1:km
       e2(n,1)=(2*n)^2;
   end
   s3=spdiags([e1,e2,e3],[-1,0,1],n,n);
 
  [v,d]=eig(full(s3));
 fc= v(:,m/2);
  
    end
if real(fc(1))<0
      fc=-fc;
end

return;
end