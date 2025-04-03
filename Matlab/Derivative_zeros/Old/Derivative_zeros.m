clc; clear; close all;

% Define parameters
v = 10;
AL = 0; AR = 20;
BL = 0; BR = AR;

% Define Mathieu cosine function (Fourier Series Approximation)
mathieu_ce = @(n, q, z) compute_mathieu_ce(n, q, z);

% Define Mathieu sine function (Fourier Series Approximation)
mathieu_se = @(n, q, z) compute_mathieu_se(n, q, z);

% Compute numerical derivatives using central finite differences
mathieu_ce_prime = @(n, q, z) (mathieu_ce(n, q, z + 1e-5) - mathieu_ce(n, q, z - 1e-5)) / (2e-5);
mathieu_se_prime = @(n, q, z) (mathieu_se(n, q, z + 1e-5) - mathieu_se(n, q, z - 1e-5)) / (2e-5);

% Define F and G functions
F = @(n, A, B) 1i * mathieu_ce_prime(n, 1/8 * (A.^2 - B.^2), 1i * atanh(B./A));
G = @(n, A, B) mathieu_se_prime(n, 1/8 * (A.^2 - B.^2), 1i * atanh(B./A));

% Define grid
A_vals = linspace(AL, AR, 50);
B_vals = linspace(BL, BR, 50);
[A, B] = meshgrid(A_vals, B_vals);


% Compute function values
F1_vals = real(F(1, A, B));
G1_vals = real(G(1, A, B));
F2_vals = real(F(2, A, B));
G2_vals = real(G(2, A, B));

% Create figures
figure;


subplot(2,2,1);
plotter(A,B,F1_vals)

subplot(2,2,2);
plotter(A,B,G1_vals)

subplot(2,2,3);
plotter(A,B,F2_vals)

subplot(2,2,4);
plotter(A,B,G2_vals)


% Function Definitions

function ce_val = compute_mathieu_ce(n, q, x)
    % Computes the derivative of cen(x, q, n) using finite differences
    h = 1e-5;
    ce_val = (cen(x + h, q, n) - cen(x - h, q, n)) / (2 * h);
end

function se_val = compute_mathieu_se(n, q, x)
    % Computes the derivative of sen(x, q, n) using finite differences
    h = 1e-5;
    se_val = (sen(x + h, q, n) - sen(x - h, q, n)) / (2 * h);
end



function plotter(A,B,F)
surf(A, B, F, 'EdgeColor', 'none'); hold on;
contour3(A, B, F, [0 0], 'k', 'LineWidth', 2);
axis tight;
caxis([-10 10])
view([0 90])
end