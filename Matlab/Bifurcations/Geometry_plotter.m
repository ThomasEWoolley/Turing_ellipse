ccc

p=[]; 
p=stanparam(p);

pde=EllipsePDE(2,1,10,100,0.9);
p.np=pde.grid.nPoints;  p.pdeo=pde; 
u=2*ones(p.np,1);
p.u=u;

plotsol(p,1,1,1)
view(2)
axis equal