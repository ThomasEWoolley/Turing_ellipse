%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   extract_one_value FUNCTION     
%   value=extract_one_value(KF,t,vncoeffs)
%
%   INPUTS:     -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd
%                
%               -t= values of order n:   
%                          KF=1   t=[0 2 4 6 ... (2*ncoeffs-2)]
%                          KF=2   t=[1 3 5 7 ... (2*ncoeffs-1)]
%                          KF=3   t=[2 4 6 8 ... (2*ncoeffs)]
%                          KF=4   t=[1 3 5 7 ... (2*ncoeffs-1)]
%               -vncoeffs= vector computed for all ncoeffs values
%
%   OUTPUT:     -value= extracted single value, corresponding to t 
%                     
%   Default value for number of coefficients: ncoeffs=25 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   extract_one_value FUNCTION CALL
function value=extract_one_value(KF,t,vncoeffs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-even---KF == 1----------------------------------------

if KF == 1  
%    ik = 0:24;
%    vt=2*ik;        % even values of order n
    in=fix(t/2)+1;
    value=vncoeffs(in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------even-odd---KF == 2 ----------------------------------------

elseif KF == 2  
%    ik = 0:24;
%    vt=2*ik+1;           % odd values of order n
    in=fix((t-1)/2)+1;
    value=vncoeffs(in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------odd-even--KF == 3 -----------------------------------------
elseif KF == 3    
%    ik=1:25; 
%    vt=2*ik;     % even values of order n
    in=fix(t/2);
    value=vncoeffs(in);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %---------------odd-odd--KF == 4 -----------------------------------------
 elseif KF == 4  
%    ik=0:24;
%    vt=2*ik+1;    % odd values of order n
    in=fix((t-1)/2)+1;
    value=vncoeffs(in);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 