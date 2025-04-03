clc; clear; close all;

% Define parameters
v = 10;
AL = 0; AR = 30;
BL = 0; BR = AR;

% Define the function to plot (Modify this only!)
F = @(n, A, B) 1i*G_dmathieuq(1,n,1/8*(A^2-B^2),1i*atanh(B./A));
%mathieu_ce_prime(n, 1/8*(A^2-B^2), 1i*atanh(B./A));
G = @(n, A, B) G_dmathieuq(2,n,1/8*(A^2-B^2),1i*atanh(B./A));


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



function plotter(A,B,F)
    % Generic plotting function for surf & contour
    I=((A)>B);
    F(~I)=nan;
    F(B==0)=nan;
    pcolor(A.*I, B.*I, F.*I); shading interp;
    hold on;
    I=((A-0.5)>B);
    F(~I)=nan;
    contour(A.*I, B.*I, F.*I, [0 0], 'k', 'LineWidth', 2);
    axis tight;
    caxis([-10 10])
    view([0 90])
    xlabel('$A$')
    ylabel('$B$')
end
