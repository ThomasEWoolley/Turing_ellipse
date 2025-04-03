%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DERIVATIVE OF RADIAL MATHIEU FUNCTION OF THE FIRST KIND  
%
%   y = dJpm(KF,u,q,mv,nmax)    [p,m = e,o (even,odd)]
%   
%   INPUTS:     -u= value of radial coordinate to compute function 
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients
%               -nmax= maximum order 
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd                           
%   OUTPUTS:    -y= vector of derivative values for all 'nmax' orders
%                   
%   'mv' is determined beforehand with function 'eig_Spm'
%   The Radial Mathieu Function is approximated by an expansion
%   of product of Bessel functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   dJpm FUNCTION CALL
function y = dJpm(KF,u,q,mv,nmax)

    v1=sqrt(q)*exp(-u);   
    v2=sqrt(q)*exp(u);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1  
     
    ik=0:24;  vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A0=Apm(1);
      
    yc = A0*(v1*besselj(1,v1)*besselj(0,v2)-v2*besselj(0,v1)*besselj(1,v2));
     
            for j=2:ncoeffs
             
    jm=fix(j-1);
    yc = yc + ((-1)^jm)*Apm(j)*(v1*besselj(j,v1)*besselj(jm,v2)- ...
                      v2*besselj(jm,v1)*besselj(j,v2));
           end         

    r=fix(in/2);
    coef=((-1)^r)*sqrt(pi/2)/A0;
    yc=yc*coef;
    
    y=[y yc];
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd---KF == 2-----------------------------------------

elseif KF == 2  
    
    ik=0:24;  vt=2*ik+1;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A1=Apm(1);
          
    yc = A1*((v2-v1)*(besselj(0,v1)*besselj(0,v2)-besselj(1,v2)*besselj(1,v1))+...
    besselj(1,v1)*besselj(0,v2)-besselj(0,v1)*besselj(1,v2));
     
            for j=2:ncoeffs
             
    jm=fix(j-1);
    jm1=fix(2*j-1);  
    
    yc = yc + ((-1)^jm)*Apm(j)*((v2-v1)*(besselj(jm,v1)*besselj(jm,v2)- ...
      besselj(j,v1)*besselj(j,v2))+jm1*(besselj(j,v1)*besselj(jm,v2)- ...
             besselj(jm,v1)*besselj(j,v2)));
         
            end
     r=fix((in-1)/2);
     coef=((-1)^r)*sqrt(pi/2)/A1;
     yc=yc*coef;
     
     y=[y yc];
     
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-even--KF == 3------------------------------------------

elseif KF == 3
    
    ik=1:25;  vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A2=Apm(1);
    
    yc = -A2*4*(besselj(0,v1)*besselj(0,v2)+ ...
   (cosh(2*u))*besselj(1,v1)*besselj(1,v2)-(1/v1)*besselj(1,v1)* ...
   besselj(0,v2)-(1/v2)*besselj(0,v1)*besselj(1,v2));
     
            for j=2:ncoeffs
             
    jm=fix(j-1);
    
    yc = yc + ((-1)^j)*Apm(j)*4*j*(besselj(jm,v1)*besselj(jm,v2)+ ...
       (cosh(2*u))*besselj(j,v1)*besselj(j,v2)-j*((1/v1)*besselj(j,v1)* ...
       besselj(jm,v2)+(1/v2)*besselj(jm,v1)*besselj(j,v2)));  

           end

   r=fix(in/2);
   coef=((-1)^r)*sqrt(pi/2)/A2;
   yc=yc*coef;
  
   y=[y yc];
   end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-odd--KF == 4-------------------------------------------

elseif KF == 4
   
    ik=0:24;    vt=2*ik+1;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);  in=vt(k);
    A1=Apm(1);
    
    yc=A1*((v1+v2)*(besselj(0,v1)*besselj(0,v2)+besselj(1,v2)*besselj(1,v1))-...
                    besselj(1,v1)*besselj(0,v2)-besselj(0,v1)*besselj(1,v2));
     
            for j=2:ncoeffs
             
    jm=fix(j-1);
    jm1=fix(2*j-1);   
    
    yc = yc + ((-1)^jm)*Apm(j)*((v1+v2)*(besselj(jm,v1)*besselj(jm,v2)+ ...
      besselj(j,v2)*besselj(j,v1))-jm1*(besselj(j,v1)*besselj(jm,v2)+ ...
      besselj(jm,v1)*besselj(j,v2)));
  
            end

    r=fix((in-1)/2);
    coef=((-1)^r)*sqrt(pi/2)/A1;
    yc=yc*coef;
    
    y=[y yc];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end