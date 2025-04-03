function [y_deriv, a_val] = mathieu_ce_deriv(n, q, z, N)
% mathieu_ce_deriv computes the derivative d/dz of the Mathieu cosine function
% MathieuCE(n,q,z) for any nonnegative integer order n.
%
% In Maple’s convention the function is normalized so that:
%   MathieuCE(n,q,0) = 1.
%
% The Fourier–series representations used here are:
%
% For even orders (n even):
%   ce(z) = sum_{m=0}^{N} A_m*cos(2*m*z),
% with the recurrence:
%   (a - 2q)A_0 - 2q*A_1 = 0,
%   -q*A_{m-1} + (a - 4*m^2)A_m - q*A_{m+1} = 0,  for m >= 1.
%
% Its derivative is computed by:
%   ce'(z) = - sum_{m=0}^{N} (2*m)*A_m*sin(2*m*z),
% and the coefficients A_m are normalized so that sum(A) = ce(0) = 1.
%
% For odd orders (n odd):
%   ce(z) = sum_{m=0}^{N} B_m*cos((2*m+1)*z),
% with the recurrence:
%   (a - 1)B_0 - q*B_1 = 0,
%   -q*B_{m-1} + (a - (2*m+1)^2)B_m - q*B_{m+1} = 0,  for m >= 1.
%
% Its derivative is:
%   ce'(z) = - sum_{m=0}^{N} (2*m+1)*B_m*sin((2*m+1)*z),
% with the normalization sum(B) = ce(0) = 1.
%
% Input:
%   n   - nonnegative integer order.
%   q   - Mathieu parameter.
%   z   - evaluation point(s) in radians (scalar or vector).
%   N   - truncation order (default is 10; use a higher N for improved accuracy).
%
% Output:
%   y_deriv - the derivative d/dz MathieuCE(n,q,z), normalized so that ce(0)=1.
%   a_val   - the characteristic value (eigenvalue) corresponding to the chosen mode.
%
% Example:
%   [yp, a] = mathieu_ce_deriv(1, 3, 2, 200);
%   % For n = 1, q = 3, and z = 2, Maple gives
%   % evalf(MathieuCPrime(1,3,2)) ≈ -0.99483.
%
% Author: [Your Name]
% Date: [Today's Date]

if nargin < 4
    N = 10;
end

if n < 0 || floor(n) ~= n
    error('Order n must be a nonnegative integer.');
end

if mod(n,2) == 0
    %% EVEN ORDER: n = 0, 2, 4, ...
    % Construct the (N+1)x(N+1) matrix for even Mathieu functions.
    M = zeros(N+1, N+1);
    % Row for m = 0: (a - 2q)A0 - 2q*A1 = 0  <=>  a*A0 = 2q*(A0 + A1)
    M(1,1) = 2*q;
    if N >= 1
        M(1,2) = 2*q;
    end
    for m = 1:(N-1)
        i = m + 1;
        M(i, i-1) = q;
        M(i, i)   = 4*m^2;
        M(i, i+1) = q;
    end
    if N >= 1
        i = N + 1;
        M(i, i-1) = q;
        M(i, i)   = 4*N^2;
    end
    % Solve the eigenvalue problem M*A = a*A.
    [V, D] = eig(M);
    eigvals = diag(D);
    [eigvals_sorted, ind] = sort(eigvals);
    V = V(:, ind);
    % For even Mathieu functions, choose mode_index = n/2 + 1.
    mode_index = n/2 + 1;
    if mode_index > size(V,2)
        error('Order n is too large for the given truncation order N.');
    end
    a_val = eigvals_sorted(mode_index);
    A = V(:, mode_index);
    % Normalize A so that sum(A) = 1 (i.e. ce(0) = 1).
    S = sum(A);
    A = A / S;
    % Compute the derivative: ce'(z) = - sum_{m=0}^{N} (2*m)*A(m+1)*sin(2*m*z).
    y = zeros(size(z));
    for m = 0:N
        y = y - (2*m)*A(m+1)*sin(2*m*z);
    end
    y_deriv = y;
    
else
    %% ODD ORDER: n = 1, 3, 5, ...
    % Construct the (N+1)x(N+1) matrix for odd Mathieu functions.
    % The recurrence is: (a - 1)*B0 - q*B1 = 0, and for m>=1:
    % -q*B_{m-1} + (a - (2*m+1)^2)*B_m - q*B_{m+1} = 0.
    M = zeros(N+1, N+1);
    % Row for m = 0:
    M(1,1) = 1;
    if N >= 1
        M(1,2) = q;
    end
    for m = 1:(N-1)
        i = m + 1;
        M(i, i-1) = q;
        M(i, i)   = (2*m+1)^2;
        M(i, i+1) = q;
    end
    if N >= 1
        i = N + 1;
        M(i, i-1) = q;
        M(i, i)   = (2*N+1)^2;
    end
    [V, D] = eig(M);
    eigvals = diag(D);
    [eigvals_sorted, ind] = sort(eigvals);
    V = V(:, ind);
    % For odd Mathieu functions, choose mode_index = (n+1)/2.
    mode_index = (n+1)/2;
    if mode_index > size(V,2)
        error('Order n is too large for the given truncation order N.');
    end
    a_val = eigvals_sorted(mode_index);
    B = V(:, mode_index);
    % Normalize B so that sum(B) = 1 (i.e. ce(0) = 1).
    S = sum(B);
    B = B / S;
    % Compute the derivative:
    % ce(z) = sum_{m=0}^{N} B_m*cos((2*m+1)*z)  so that
    % ce'(z) = - sum_{m=0}^{N} (2*m+1)*B_m*sin((2*m+1)*z).
    y = zeros(size(z));
    for m = 0:N
        y = y - (2*m+1)*B(m+1)*sin((2*m+1)*z);
    end
    y_deriv = y;
    
end

end
