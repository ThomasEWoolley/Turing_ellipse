function [y, a_val] = mathieu_ce(mode, q, theta, N)
% mathieu_ce calculates the Mathieu cosine function ce for a given mode,
% parameter q and angle theta.
%
% This function computes the even Mathieu function defined by the series
% expansion:
%
%   ce(theta) = sum_{m=0}^{N} A_m * cos(2*m*theta),
%
% where the coefficients A_m satisfy the recurrence:
%   a*A_0 - 2*q*A_1 = 0,               (m = 0)
%   -q*A_{m-1} + (a - 4*m^2)*A_m - q*A_{m+1} = 0,  (m >= 1).
%
% The characteristic value a and the coefficient vector A are determined by
% solving the eigenvalue problem M*A = a*A, where M is a (N+1)x(N+1)
% tridiagonal matrix.
%
% Input:
%   mode  - the mode number (0 gives the lowest even mode, 1 the next, etc.)
%   q     - the Mathieu parameter.
%   theta - the angle (in radians) at which to evaluate the function.
%           This can be a scalar or a vector.
%   N     - truncation number (optional, default is 10).
%
% Output:
%   y     - the value of the Mathieu cosine function at theta.
%   a_val - the characteristic value (eigenvalue) corresponding to the mode.
%
% Example:
%   [y, a] = mathieu_ce(0, 1, linspace(0, 2*pi, 100), 20);
%   plot(linspace(0, 2*pi, 100), y);

if nargin < 4
    N = 10; % default truncation number
end

% Construct the tridiagonal matrix M of size (N+1)x(N+1).
% The series expansion uses indices m = 0,1,...,N.
% For m = 0, the recurrence is: a*A_0 - 2*q*A_1 = 0.
% For m >= 1:  -q*A_{m-1} + (a - 4*m^2)*A_m - q*A_{m+1} = 0.
M = zeros(N+1, N+1);

% For m = 0 (index i = 1)
M(1,1) = 0;
if N >= 1
    M(1,2) = -2*q;
end

% For m = 1 to N-1 (indices i = 2 to N)
for m = 1:(N-1)
    i = m + 1;
    M(i, i-1) = -q;
    M(i, i)   = 4*m^2;
    M(i, i+1) = -q;
end

% For m = N (index i = N+1)
if N >= 1
    i = N + 1;
    M(i, i-1) = -q;
    M(i, i)   = 4*N^2;
end

% Solve the eigenvalue problem.
[V, D] = eig(M);
eigvals = diag(D);

% Sort eigenvalues (and corresponding eigenvectors) in ascending order.
[eigvals_sorted, ind] = sort(eigvals);
V = V(:, ind);

% Check that the mode index is valid.
if mode < 0 || mode > N
    error('Mode must be between 0 and N (with N the truncation order).');
end

% For the even Mathieu function, mode 0 corresponds to the smallest eigenvalue.
selected_index = mode + 1;
a_val = eigvals_sorted(selected_index);
A = V(:, selected_index);

% Normalise the coefficient vector.
if abs(A(1)) > 1e-12
    A = A / A(1);
else
    A = A / norm(A);
end

% Compute ce(theta) = sum_{m=0}^{N} A(m+1)*cos(2*m*theta).
y = zeros(size(theta));
for m = 0:N
    y = y + A(m+1) * cos(2*m*theta);
end

end
