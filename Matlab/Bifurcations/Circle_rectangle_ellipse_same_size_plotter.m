%%
ccc
figure('Units','normalized','Position',[1 0 1 1/3])

Ls=[220 255 265 270];
l=1;
for i=Ls
    subplot(1,5,l)
M=mphload(['../../Comsol/Circle_vs_rectangle_all_same_',num2str(i),'.mph'],'Model');
Q=mphsolutioninfo(M);

hold on
% figure('Units','Normalized','position',[0 0.05 1/3 0.85])

U=mpheval(M,'u','t',10000);
h=trisurf(triangulation(double(U.t+1)',U.p(1,:)',U.p(2,:)',U.d1(1,:)'));
axis equal
% axis off
view(2)
title(['$A=$',num2str(i/100)])
shading interp
if l==1
    caxis([0.999 1.001])
end
h = gca;
h.YAxis.Visible = 'off';
xlim([-3 3])
 xticks(-3:1:3)
 xtickangle(0)
 
colorbar
l=l+1;
set(gca,'fontsize',15)
end


export_fig(['../../Pictures/Circle_rectangle_ellipse_same_size.png'],'-r300')