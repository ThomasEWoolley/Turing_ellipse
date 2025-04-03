ccc
A=5;
B=1;
a=sqrt(A^2-B^2);
v=linspace(0,2*pi);
uc=atanh(B/A);
u=linspace(0,uc,2);

[U,V]=meshgrid(uc,v);
x=a*cosh(U).*cos(V);
y=a*sinh(U).*sin(V);

scatter(2*x,y,'b','filled')
hold on
A=10;
B=1;
a=sqrt(A^2-B^2);
v=linspace(0,2*pi);
uc=atanh(B/A);
u=linspace(0,uc,2);

[U,V]=meshgrid(uc,v);
x=a*cosh(U).*cos(V);
y=a*sinh(U).*sin(V);

scatter(x,y,'r')