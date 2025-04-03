function  w = dbesselj(nu,z,scale)
%     ==============================================================
%     Purpose: Compute the derivative of Bessel function 
%     
%     Input:   nu--- Order of Bessel functions
%     z  ---Argument of Bessel functions
%     Output:   w--- the derivative of Bessel function
%
%    Editor: Da Ma, Southeast University, China.
%     ==============================================================
if nargin == 2, scale = 0; end
if nu==0
    w=-besselj(1,z)
else
    w=0.5.*(besselj(nu-1,z)-besselj(nu+1,z));
end


