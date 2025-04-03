function csd = G_dmathieuq4(kf, m, q, xr, varargin)
% G_DMATHIEUQ2  Compute the derivative of angular Mathieu functions with improved accuracy.
%
%   csd = G_dmathieuq2(kf, m, q, xr)
%   csd = G_dmathieuq2(kf, m, q, xr, 'MaxTerms', N, 'Tol', tol)
%
%   Inputs:
%       kf    - Function code:
%               1 for even angular Mathieu function derivative (cem')
%               2 for odd angular Mathieu function derivative (sem')
%       m     - Order of the Mathieu function.
%       q     - Parameter of the Mathieu function (can be positive, negative,
%               or complex). Can be a vector.
%       xr    - Argument of the Mathieu function (in radians).
%
%   Optional Name-Value Pairs:
%       'MaxTerms' - Maximum number of terms in the series expansion 
%                    (default: auto-determined from q and m).
%       'Tol'      - Tolerance for series truncation (default: 1e-10).
%
%   Output:
%       csd   - The derivative of the Mathieu function evaluated at xr.
%
%   This version has been modified to use a more robust eigenvector normalization
%   (forcing ce(0)=1) and to adapt the number of series terms for large |q|.
%
%   Author: [Your Name]
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
    % Use an adaptive heuristic for maximum number of terms:
    % For small |q| use a modest number, for large |q| add extra terms.
    maxCol = ceil(20 + 10*sqrt(max(abs(q))));
else
    maxCol = maxTerms;
end
% Preallocate coefficient array (each row corresponds to one q value).
fg = zeros(lq, maxCol);

% Determine kd based on kf and parity of m.
if (kf == 1 && mod(m,2) == 0)
    kd = 1;  % even function derivative (cem') for even m
elseif (kf == 1 && mod(m,2) ~= 0)
    kd = 2;  % even function derivative (cem') for odd m
elseif (kf == 2 && mod(m,2) ~= 0)
    kd = 3;  % odd function derivative (sem') for odd m
elseif (kf == 2 && mod(m,2) == 0)
    kd = 4;  % odd function derivative (sem') for even m
else
    error('Invalid combination of kf and m.');
end

% Determine series truncation index (km) for each q.
qm = zeros(lq,1);
km = zeros(lq,1);
for k = 1:lq
    if abs(q(k)) <= 1
        % Heuristic for small |q|
        qm(k) = 7.5 + 56.1*sqrt(abs(q(k))) - 134.7*abs(q(k)) + 90.7*sqrt(abs(q(k)))*abs(q(k));
    else
        % For larger |q| use a different heuristic and add extra terms.
        qm(k) = 17.0 + 3.1*sqrt(abs(q(k))) - 0.126*abs(q(k)) + 0.0037*sqrt(abs(q(k)))*abs(q(k));
    end
    km(k) = fix(qm(k) + 0.5*fix(m));
    % For very large |q|, add extra terms.
    if abs(q(k)) > 10
        km(k) = km(k) + 10;
    end
    if abs(q(k)) >= 30
        km(k) = km(k) + 20;  % extra terms for very large q
    end
    % Respect a maximum if provided.
    if ~isempty(maxTerms) && km(k) > maxTerms
        km(k) = maxTerms;
    end
end

% Compute expansion coefficients for each q using the improved routine.
for k = 1:lq
    fc = fcoef1(kd, m, q(k), km(k));
    fc = fc(:).';  % ensure row vector
    numCoeffs = length(fc);
    fg(k, 1:numCoeffs) = fc;
end

% Determine starting index for possible early termination.
ic = fix(m/2) + 1;

% Disable early termination if parameters are large.
disableBreak = (max(abs(q)) >= 30) || (abs(xr) >= 30);

% Compute the derivative via series summation using Kahan summation.
csd = zeros(1, lq);
for vm = 1:lq
    if isreal(q(vm)) && real(q(vm)) < 0
        csd(vm) = NaN;
    else
        sumVal = 0;
        comp   = 0; % compensation for lost low-order bits
        % Sum terms from k = 1 to km(vm)
        for k = 1:km(vm)
            switch kd
                case 1  % even function, even order: series in cos(2*(k-1)*x)
                    freq = 2*(k-1);
                    term = - freq * fg(vm, k) * sin(freq*xr);
                case 2  % even function, odd order: series in cos((2*(k-1)+1)*x)
                    freq = 2*(k-1)+1;
                    term = - freq * fg(vm, k) * sin(freq*xr);
                case 3  % odd function, odd order: series in cos((2*(k-1)+1)*x)
                    freq = 2*(k-1)+1;
                    term =  freq * fg(vm, k) * cos(freq*xr);
                case 4  % odd function, even order: series in cos(2*k*x)
                    freq = 2*k;
                    term =  freq * fg(vm, k) * cos(freq*xr);
            end
            % Kahan summation update
            y = term - comp;
            t = sumVal + y;
            comp = (t - sumVal) - y;
            sumVal = t;
            
            % For low-parameter values, allow early termination.
            if ~disableBreak && k >= ic && abs(term) < tol * abs(sumVal)
                break;
            end
        end
        csd(vm) = sumVal;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc = fcoef1(kd, m, q, km)
% FCOEF1  Compute expansion coefficients for angular Mathieu functions.
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
%       km - Number of terms in the series expansion.
%
%   Output:
%       fc - Expansion coefficients (row vector) normalized so that the sum equals 1.
%
%   This version uses a full (dense) eigenvalue solver and normalizes the
%   eigenvector so that ce(0)=sum(coefficients)=1.
%
%   Author: [Your Name]
%   Date: [Date]

n = km;
switch kd
    case 1
        % Even Mathieu function (cem) for even m.
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*(j-1))^2;
        end
        % Build full tridiagonal matrix.
        A = full(spdiags([e1, e2, e3], [-1, 0, 1], n, n));
        % Adjust first two rows to impose correct recurrence.
        A(1,1) = 1e-14;  % small number to avoid singularity
        if n >= 2
            A(2,1) = 2*q;
        end
    case 2
        % Even Mathieu function (cem) for odd m.
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*j - 1)^2;
        end
        A = full(spdiags([e1, e2, e3], [-1, 0, 1], n, n));
        A(1,1) = 1 + q;
    case 3
        % Odd Mathieu function (sem) for odd m.
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*j - 1)^2;
        end
        A = full(spdiags([e1, e2, e3], [-1, 0, 1], n, n));
        A(1,1) = 1 - q;
    case 4
        % Odd Mathieu function (sem) for even m.
        e1 = q * ones(n,1);
        e2 = zeros(n,1);
        e3 = q * ones(n,1);
        for j = 1:n
            e2(j) = (2*j)^2;
        end
        A = full(spdiags([e1, e2, e3], [-1, 0, 1], n, n));
    otherwise
        error('Invalid kd value.');
end

% Solve the eigenvalue problem.
[V, D] = eig(A);
[dSorted, sortIdx] = sort(diag(D));
V = V(:, sortIdx);

% Select the eigenvector index based on m.
switch kd
    case 1
        idx = (m + 2) / 2;
    case {2, 3}
        idx = (m + 1) / 2;
    case 4
        idx = m / 2;
end
if idx < 1 || idx > size(V,2)
    error('Eigenvector index out of range. Increase km.');
end
fc = V(:, idx);
% Normalize the coefficients so that their sum is 1.
normFactor = sum(fc);
if normFactor == 0
    normFactor = 1;
end
fc = fc / normFactor;
% Ensure the first coefficient is positive.
if real(fc(1)) < 0
    fc = -fc;
end

end
