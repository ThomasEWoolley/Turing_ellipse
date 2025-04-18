%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DERIVATIVE OF ANGULAR MATHIEU FUNCTION 
%   y = dSpm(KF,v,mv,nmax)     [p,m = e,o (even,odd)]
%   
%   INPUTS:     -v= value of angular coordinate in radians 
%               -mv= matrix of expansion coefficients
%               -nmax= maximum order
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd                                                        
%   OUTPUTS:    -y= vector of derivative values for all 'nmax' orders
%   'mv' is determined beforehand with function 'eig_Spm'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   dSpm FUNCTION CALL
function y = dSpm(KF,v,mv,nmax)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1-------------------------------------

if KF == 1  

    ik=0:24;   vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k); 
   
     yc = 0;
            for j = 2:ncoeffs
                jc = fix(2*(j-1));
                yc = yc - jc*Apm(j)*sin(jc*v);
            end
            
     y=[y yc];
     
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd---KF == 2-----------------------------------------

elseif KF == 2  

    ik=0:24;  vt=2*ik+1;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k); 
   
    yc = 0;
            for j = 1:ncoeffs
                yc = yc - (2*j-1)*Apm(j)*sin((2*j-1)*v);
            end
    
     y=[y yc];
     
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-even---KF == 3-----------------------------------------

elseif KF == 3  

    ik=1:25;  vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k); 
   
    yc = 0;
    
            for j = 1:ncoeffs               
                yc = yc + 2*j*Apm(j)*cos(2*j*v);
            end
    
     y=[y yc]; 
     
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-odd---KF == 4------------------------------------------

elseif KF == 4  

    ik=0:24;  vt=2*ik+1; 
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k); 
   
    yc = 0;
            for j = 1:ncoeffs             
                yc = yc + (2*j-1)*Apm(j)*cos((2*j-1)*v);
            end
     
     y=[y yc]; 
     
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

