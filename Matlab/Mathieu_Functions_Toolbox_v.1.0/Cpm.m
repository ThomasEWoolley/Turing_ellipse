%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CORRELATION FACTOR    
%   y = Cpm(KF,mv,mvv,nmax)   [p,m = e,o (even,odd)]
%
%   INPUTS:     -mv= matrix of expansion coefficients
%               -mvv= matrix of expansion coefficients
%               -nmax= maximum order
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd                                                        
%   OUTPUTS:    -y= vector of correlation factors for all 'nmax' orders                   
%   'mv', 'mvv' are determined beforehand with function 'eig_Spm'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Cpm FUNCTION CALL
function y = Cpm(KF,mv,mvv,nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1  

    ik=0:24;  vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
    
    Apm=mv(:,k);   AApm=mvv(:,k);
    A0=Apm(1);     AA0=AApm(1);
    
         yc = 2*pi*A0*AA0;
            for j = 2:ncoeffs
                yc = yc + pi * Apm(j) * AApm(j);
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
    
    Apm=mv(:,k);  AApm=mvv(:,k);    
    
    yc = pi * sum(Apm.*AApm);

    y=[y yc];
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------odd-even--KF == 3-----------------------------------------

elseif KF == 3  
    
    ik=1:25;  vt=2*ik;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
           
    Apm=mv(:,k);  AApm=mvv(:,k);     
    
    yc = pi * sum(Apm.*AApm);

    y=[y yc];
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------odd-odd--KF == 4------------------------------------------

elseif KF == 4  
    
    ik=0:24;  vt=2*ik+1;
    ncoeffs=length(vt);
    
    y=[];
    
    for k=1:nmax
           
    Apm=mv(:,k);  AApm=mvv(:,k);     
    
    yc = pi * sum(Apm.*AApm);

    y=[y yc];
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
    