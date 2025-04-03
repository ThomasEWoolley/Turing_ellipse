ccc

cols=[1 0 0;
      0 0.5 0;
      0 0 1;
      200/256, 200/256, 0
      203/256, 195/256, 227/256
      0 0 0];

B=1;

A0=2.66;
F = @(A, B) 1i*G_dmathieuq(1,1,1/8*(A^2-B^2),1i*atanh(B./A));
A=fsolve(@(A)F(A,B),A0,optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt'));
a=sqrt(A^2-B^2);
uc=atanh(B/A);

[u,v]=meshgrid(linspace(0,uc),linspace(0,2*pi));

z=arrayfun(@(u)G_mathieuq(1,1,1/8*(A^2-B^2),u*1i),u).*arrayfun(@(v)G_mathieuq(1,1,1/8*(A^2-B^2),v),v);
pcolor(a*cosh(u).*cos(v),a*sinh(u).*sin(v),z)
shading interp

axis equal
axis(5*[-1 1 -1 1])
colorbar
xlabel('$x$')
ylabel('$y$')
set(gca,'fontsize',20)
export_fig('../../Pictures/Ce','-r300')
Ac=A;
Bc=B;
%%
figure
A=4;
B0=2.8;
F = @(A, B) G_dmathieuq(2,1,1/8*(A^2-B^2),1i*atanh(B./A));
B=fsolve(@(B)F(A,B),B0,optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt'));
a=sqrt(A^2-B^2);
uc=atanh(B/A);

[u,v]=meshgrid(linspace(0,uc),linspace(0,2*pi));

z=arrayfun(@(u)-1i*G_mathieuq(2,1,1/8*(A^2-B^2),1i*u),u).*arrayfun(@(v)G_mathieuq(2,1,1/8*(A^2-B^2),v),v);
pcolor(a*cosh(u).*cos(v),a*sinh(u).*sin(v),z)
shading interp
axis equal
axis(5*[-1 1 -1 1])
colorbar
xlabel('$x$')
ylabel('$y$')
set(gca,'fontsize',20)
export_fig('../../Pictures/Se','-r300')
As=A;
Bs=B;
save('Ce_Se.mat','As','Bs','Ac',"Bc")

