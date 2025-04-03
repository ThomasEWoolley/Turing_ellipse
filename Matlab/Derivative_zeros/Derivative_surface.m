clc; clear; close all;

% Define parameters
v = 10;
AL = 0; AR = 15;
BL = 0; BR = AR;

% Define the function to plot (Modify this only!)
F = @(n, A, B) 1i*G_dmathieuq2(1,n,1/8*(A^2-B^2),1i*atanh(B./A));
%mathieu_ce_prime(n, 1/8*(A^2-B^2), 1i*atanh(B./A));
G = @(n, A, B) G_dmathieuq2(2,n,1/8*(A^2-B^2),1i*atanh(B./A));


% Define grid
nn=200;
A_vals = linspace(AL, AR, nn);
B_vals = linspace(BL, BR, nn);
[A, B] = meshgrid(A_vals, B_vals);

% Compute function values (use arrayfun to apply element-wise)
F1_vals = real(arrayfun(@(a, b) F(1, a, b), A, B));
G1_vals = real(arrayfun(@(a, b) G(1, a, b), A, B));
F2_vals = real(arrayfun(@(a, b) F(2, a, b), A, B));
G2_vals = real(arrayfun(@(a, b) G(2, a, b), A, B));

%%
close all
% Create figures
load Ce_Se.mat
figure;

subplot(2,2,1);
plotter(A,B,F1_vals)
title('$\textrm{d} Ce_1/\textrm{d} u$')
plot(Ac,Bc,'ko')
subplot(2,2,2);
plotter(A,B,G1_vals)
title('$\textrm{d} Se_1/\textrm{d} u$')
plot(As,Bs,'ko')

subplot(2,2,3);
plotter(A,B,F2_vals)
title('$\textrm{d} Ce_2/\textrm{d} u$')

subplot(2,2,4);
plotter(A,B,G2_vals)
title('$\textrm{d} Se_2/\textrm{d} u$')


export_fig('../../Pictures/Derivative_zeros','-r300')

function plotter(A,B,F)
    % Generic plotting function for surf & contour
    I=((A)>B);
    F(~I)=nan;
    F(B==0)=nan;
    pcolor(A.*I, B.*I, F.*I); shading interp;
    hold on;
    I=((A-0.2)>B);
    F(~I)=nan;
    contour(A.*I, B.*I, F.*I, [0 0], 'k', 'LineWidth', 2);
    axis tight;
    caxis([-1 1])
    view([0 90])
    xlabel('$A$')
    ylabel('$B$')
    set(gca,'fontsize',12)
end
