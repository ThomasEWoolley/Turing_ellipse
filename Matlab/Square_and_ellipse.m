ccc
theta=linspace(0,2*pi);
r=2.6;
plot(r*cos(theta),r*sin(theta),LineWidth=3)
hold on
plot(pi/sqrt(2)*[-1 1 1 -1 -1],pi/sqrt(2)*[1 1 -1 -1 1],LineWidth=3)
axis equal