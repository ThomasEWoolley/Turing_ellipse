ccc
cols=[1 0 0;
      0 0.5 0;
      0 0 1;
      200/256, 200/256, 0
      203/256, 195/256, 227/256
      0 0 0];

B=linspace(2.6,0,1e3);
figure('Units','normalized','Position',[1 0 1 1/3])
subplot(1,5,1)
A0=2.66;
F = @(A, B) 1i*G_dmathieuq4(1,1,1/8*(A^2-B^2),1i*atanh(B./A));

i=1;
for BB=B

a(i) = fsolve(@(A)F(A,BB),A0,optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt'));
A0=a(i);
i=i+1;
end

plot(a,B,'k','linewidth',3)
hold on
plot(pi/2*sqrt(2),0,'s','linewidth',3,'MarkerFaceColor',cols(6,:),'MarkerEdgeColor',cols(6,:))
plot(2.6022,2.6022,'o','linewidth',3,'MarkerFaceColor',cols(5,:),'MarkerEdgeColor',cols(5,:))
plot(3*[0 1],3*[0,1],'--')
xlabel('$A$')
ylabel('$B$')
xlim([2.2 2.7])

Bifpoint = [ 2.6663       2.6583       2.6456       2.6289];
Points = 0.5:0.5:2;
for i=1:4
    cols(i,:)
plot(Bifpoint(i),Points(i),'d','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',10)
end
set(gca,'fontsize',12)
export_fig('../../Pictures/dce0','-r300')