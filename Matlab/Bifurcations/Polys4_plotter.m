ccc

p=[]; 
p=stanparam(p);
pde=Polys4(1,10,5,.1);
p.np=pde.grid.nPoints;  p.pdeo=pde; 

u=2*ones(p.np,1);

p.u=u;

plotsol(p,1,1,1)
view(2)
axis equal