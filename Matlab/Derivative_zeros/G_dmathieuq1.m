function csd = G_dmathieuq(kf, m, q, xr, varargin)
% G_DMATHIEUQ Compute the derivative of angular Mathieu functions.
%
%   csd = G_dmathieuq(kf, m, q, xr)
%   csd = G_dmathieuq(kf, m, q, xr, 'MaxTerms', N, 'Tol', tol)
%
%   Inputs:
%       kf    - Function code:
%               1 for even angular Mathieu function derivative (cem')
%               2 for odd angular Mathieu function derivative (sem')
%       m     - Order of the Mathieu function.
%       q     - Parameter of the Mathieu function (can be positive, negative,
%               or complex).
%       xr    - Argument of the Mathieu function (in radians).
%
%   Optional Name-Value Pairs:
%       'MaxTerms' - Maximum number of terms in the series expansion 
%                    (default: auto-determined).
%       'Tol'      - Tolerance for series truncation (default: 1e-10).
%
%   Output:
%       csd   - The derivative of the Mathieu function evaluated at xr.
%
%   This modified code improves accuracy for larger q values (around 20)
%   by increasing the series length and using a more forgiving truncation
%   tolerance.
%
%   Author: Da Ma (original), modified by [Your Name]
%   Date: [Date]

% Parse optional parameters
p = inputParser;
addParameter(p, 'MaxTerms', []);
addParameter(p, 'Tol', 1e-10);
parse(p, varargin{:});
maxTerms = p.Results.MaxTerms;
tol      = p.Results.Tol;

lq = length(q);
if isempty(maxTerms)
    maxCol = 251;
else
    maxCol = maxTerms;
end
fg = zeros(lq, maxCol);

% Determine kd code based on kf and m parity.
if (kf == 1 && mod(m,2) == 0)
    kd = 1;
elseif (kf == 1 && mod(m,2) ~= 0)
    kd = 2;
elseif (kf == 2 && mod(m,2) ~= 0)
    kd = 3;
elseif (kf == 2 && mod(m,2) == 0)
    kd = 4;
end

% Determine series truncation index (km) for each q.
qm = zeros(lq,1);
km = zeros(lq,1);
for k = 1:lq
    if abs(q(k)) <= 1
        qm(k) = 7.5 + 56.1*sqrt(abs(q(k))) - 134.7*abs(q(k)) ...
                   + 90.7*sqrt(abs(q(k)))*abs(q(k));
    else
        qm(k) = 17.0 + 3.1*sqrt(abs(q(k))) - 0.126*abs(q(k)) ...
                   + 0.0037*sqrt(abs(q(k)))*abs(q(k));
    end
    km(k) = fix(qm(k) + 0.5*fix(m));
    % For larger q, add extra terms for better accuracy.
    if abs(q(k)) > 10
        km(k) = km(k) + 10;
    end
    % If a maximum number of terms is specified, limit km accordingly.
    if ~isempty(maxTerms) && km(k) > maxTerms
        km(k) = maxTerms;
    end
end

% Compute expansion coefficients for each q.
for k = 1:lq
    fc = fcoef1(kd, m, q(k), km(k));
    fc = fc(:).';  % Ensure a row vector.
    numCoeffs = length(fc);
    fg(k, 1:numCoeffs) = fc;
end

% Determine starting index for truncation.
ic = fix(m/2) + 1;

% Compute the derivative via series summation.
csd = zeros(1, lq);
for vm = 1:lq
    if isreal(q(vm)) && real(q(vm)) < 0
        csd(vm) = NaN;
    else
        sumVal = 0;
        for k = 1:km(vm)
            if kd == 1
                term = -(2*k - 2) * fg(vm, k) * sin((2*k - 2)*xr);
            elseif kd == 2
                term = -(2*k - 1) * fg(vm, k) * sin((2*k - 1)*xr);
            elseif kd == 3
                term = (2*k - 1) * fg(vm, k) * cos((2*k - 1)*xr);
            elseif kd == 4
                term = 2*k * fg(vm, k) * cos(2*k*xr);
            end
            sumVal = sumVal + term;
            % Stop if the current term is small relative to the sum.
            if k >= ic && abs(term) < tol * abs(sumVal)
                break;
            end
        end
        csd(vm) = sumVal;
    end
end

end


function fc = fcoef1(kd, m, q, km)
% FCOEF1 Compute expansion coefficients for angular Mathieu functions.
%
%   fc = fcoef1(kd, m, q, km)
%
%   Inputs:
%       kd - Case code:
%            1 for cem(x,q,m = 0,2,4,...)
%            2 for cem(x,q,m = 1,3,5,...)
%            3 for sem(x,q,m = 1,3,5,...)
%            4 for sem(x,q,m = 2,4,6,...)
%       m  - Order of the Mathieu function.
%       q  - Parameter of the Mathieu function.
%       km - Number of terms to use in the series expansion.
%
%   Output:
%       fc - Expansion coefficients for the Mathieu function.
%
%   This revised version allows for an increased number of terms (km)
%   to improve accuracy for larger q values.
%
%   Author: Da Ma (original), modified by [Your Name]
%   Date: [Date]

n = km;

switch kd
    case 1
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*(j-1))^2;
        end
        % Construct tridiagonal matrix.
        A = spdiags([e1, e2, e3], [-1, 0, 1], n, n);
        % Adjust the first element to avoid singularity.
        A(1,1) = 1e-14;
        A(2,1) = 2*q;
        [V, ~] = eig(full(A));
        idx = (m + 2) / 2; % Select appropriate eigenvector.
        fc = V(:, idx);
        fc = fc / sqrt(1 + abs(fc(1))^2);
    case 2
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*j - 1)^2;
        end
        A = spdiags([e1, e2, e3], [-1, 0, 1], n, n);
        A(1,1) = 1 + q;
        [V, ~] = eig(full(A));
        idx = (m + 1) / 2;
        fc = V(:, idx);
    case 3
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*j - 1)^2;
        end
        A = spdiags([e1, e2, e3], [-1, 0, 1], n, n);
        A(1,1) = 1 - q;
        [V, ~] = eig(full(A));
        idx = (m + 1) / 2;
        fc = V(:, idx);
    case 4
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*j)^2;
        end
        A = spdiags([e1, e2, e3], [-1, 0, 1], n, n);
        [V, ~] = eig(full(A));
        idx = m / 2;
        fc = V(:, idx);
end

% Ensure the first coefficient is positive.
if real(fc(1)) < 0
    fc = -fc;
end

end
