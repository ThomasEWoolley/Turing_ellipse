ccc


M=mphload(['../Comsol/Circle_vs_rectangle.mph'],'Model');
Q=mphsolutioninfo(M);
%%
close all
hold on
% figure('Units','Normalized','position',[0 0.05 1/3 0.85])

U=mpheval(M,'u','t',10000);
h=trisurf(triangulation(double(U.t+1)',U.p(1,:)',U.p(2,:)',U.d1(1,:)'));
axis equal
axis off
view(2)

shading interp

caxis([0.95 1.05])
colorbar
set(gca,'fontsize',12)


export_fig(['../Pictures/Circle_rectangle_ellipse.png'],'-r300')