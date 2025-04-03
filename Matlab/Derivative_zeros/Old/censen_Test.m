clc; clear; close all;

% Define parameters
v = 10;
AL = 0; AR = 5;
BL = 0; BR = AR;

% Define the function to plot (Modify this only!)
F = @(n, A, B) mathieu_ce_prime(n, 1/8*(A^2-B^2), 1i*atanh(B./A));
G = @(n, A, B) -1i*mathieu_se_prime(n, 1/8*(A^2-B^2), 1i*atanh(B./A));

1i*F(1,5,3)
%%
% Define grid
nn=100;
A_vals = linspace(AL, AR, nn);
B_vals = linspace(BL, BR, nn);
[A, B] = meshgrid(A_vals, B_vals);

% Compute function values (use arrayfun to apply element-wise)
F1_vals = real(arrayfun(@(a, b) F(1, a, b), A, B));
G1_vals = real(arrayfun(@(a, b) G(1, a, b), A, B));
F2_vals = real(arrayfun(@(a, b) F(2, a, b), A, B));
G2_vals = real(arrayfun(@(a, b) G(2, a, b), A, B));

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


% **Function Definitions**


function ce_prime = mathieu_ce_prime(n, q, x)
    % Computes the derivative of cen(x, q, n) using finite differences
    h = 1e-5;
    ce_prime = (cen(x + h, q, n) - cen(x - h, q, n)) / (2 * h);
end

function se_prime = mathieu_se_prime(n, q, x)
    % Computes the derivative of sen(x, q, n) using finite differences
    h = 1e-5;
    se_prime = (sen(x + h, q, n) - sen(x - h, q, n)) / (2 * h);
end

function plotter(A,B,F)
    % Generic plotting function for surf & contour
    surf(A.*(A>B), B.*(A>B), F.*(A>B), 'EdgeColor', 'none'); hold on;
    contour3(A.*(A>B), B.*(A>B), F.*(A>B), [0 0], 'k', 'LineWidth', 2);
    axis tight;
    % caxis([-10 10])
    view([0 90])
    xlabel('A')
    ylabel('B')
end
