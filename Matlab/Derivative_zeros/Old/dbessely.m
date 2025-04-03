function  w = dbessely(nu,z,scale)
%     ==============================================================
%     Purpose: Compute the derivative of Bessel function of second kind
%     
%     Input:   nu--- Order of Bessel function of second kind
%     z  ---Argument of Bessel function of second kind
%     Output:   w--- the derivative of Bessel function of second kind
%
%    Editor: Da Ma, Southeast University, China.
%     ==============================================================
if nargin == 2, scale = 0; end
if nu==0
    w=-bessely(1,z)
else
    w=0.5.*(bessely(nu-1,z)-bessely(nu+1,z));
end


