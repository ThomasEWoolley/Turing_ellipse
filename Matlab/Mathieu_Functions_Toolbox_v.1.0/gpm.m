%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   JOINING FACTOR    
%   y = gpm(KF,q,mv,nmax)  [p,m = e,o (even,odd)]
%
%   INPUTS:     -q= elliptical parameter, q > 0
%               -mv= matrix of expansion coefficients
%               -nmax= maximum order
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd                                                        
%   OUTPUTS:    -y= vector of joining factors for all 'nmax' orders                   
%   'mv' is determined beforehand with function 'eig_Spm'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   gpm FUNCTION CALL
function y = gpm(KF,q,mv,nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1
    
    ik=0:24;  vt=2*ik;
%   Compute Spm at v=pi/2    
    Spm_pi2=Spm(1,pi/2,mv,nmax);                   
   
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);   in=vt(k);
    A0=Apm(1);     
    Spmk=Spm_pi2(k);                   
    r=fix(in/2);
    yc=((-1)^r)*(1/pi)*Spmk/A0; 
    y=[y yc];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd---KF == 2-----------------------------------------

elseif KF == 2  
    
    ik=0:24;  vt=2*ik+1;
%   Compute the derivative of Spm at v=pi/2
    dSpm_pi2=dSpm(2,pi/2,mv,nmax);                   
   
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);   in=vt(k);
    A1=Apm(1);     
    dSpmk=dSpm_pi2(k); 
    r=fix((in-1)/2);
    yc=-((-1)^r)*(1/(pi*sqrt(q)))*dSpmk/A1;     
    y=[y yc];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------odd-even--KF == 3-----------------------------------------

elseif KF == 3  
    
    ik=1:25;  vt=2*ik;
%   Compute the derivative of Spm at v=pi/2 
    dSpm_pi2=dSpm(3,pi/2,mv,nmax);                   
   
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);   in=vt(k);
    A2=Apm(1);     
    dSpmk=dSpm_pi2(k);
    r=fix(in/2);
    yc=((-1)^r)*(1/(q*pi))*dSpmk/A2;     
    y=[y yc];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------odd-odd--KF == 4------------------------------------------

elseif KF == 4  
    
    ik=0:24;  vt=2*ik+1;
%   Compute Spm at v=pi/2  
    Spm_pi2=Spm(4,pi/2,mv,nmax); 
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);   in=vt(k);
    A1=Apm(1);     
    Spmk=Spm_pi2(k);
    r=fix((in-1)/2);
    yc=((-1)^r)*(1/(sqrt(q)*pi))*Spmk/A1;   
    y=[y yc];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

