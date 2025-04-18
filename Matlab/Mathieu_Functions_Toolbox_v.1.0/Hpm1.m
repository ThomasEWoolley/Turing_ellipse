%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RADIAL MATHIEU FUNCTION OF THE THIRD KIND    
%   y = Hpm1(KF,u,q,mv,nmax)   [p,m = e,o (even,odd)]
%   
%   INPUTS:     -u= value of radial coordinate to compute function 
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients
%               -nmax= maximum order 
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd      
%   OUTPUTS:    -y= vector of function values for all 'nmax'orders 
%                  
%   The Radial Mathieu Function of the Third Kind is defined by analogy   
%   with the Hankel Function of the First Kind:  
%     Hpm1(KF,u,q,mv,nmax)=Jpm(KF,u,q,mv,nmax)+i*Ypm(KF,u,q,mv,nmax)
%   'mv' is determined beforehand with function 'eig_Spm'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Hpm1 FUNCTION CALL
function y = Hpm1(KF,u,q,mv,nmax)

   y1 = Jpm(KF,u,q,mv,nmax);
   y2 = Ypm(KF,u,q,mv,nmax);

   y = y1 + i*y2;

