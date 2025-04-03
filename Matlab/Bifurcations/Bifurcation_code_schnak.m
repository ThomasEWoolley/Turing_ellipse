
close all;clc;keep pphome;
B=0.5;
Fname=['Ellipsehom_',num2str(B)];

if exist(Fname, 'dir')
    rmdir([Fname,'*'],'s')
end


p=[]; p=stanparam(p); % Setting up basic p structure
p.plot.auxdict={'$L$' ,'$Du$', '$Dv$','$a$','$b$','$max u1$','$min u1$'}; % Parameter names (not needed, but nice)
par=           [ 1   , 1  ,  10  ,  1/10   , 9/10]; % Parameters used in functions and Jacobian
p.fuha.outfu=@outfn; % Output function to be plotted, names are in above.
p.nc.ilam=1; % Parameter number to continue over
p.plot.bpcmp=6; % Which component of output function to plot.
p.plot.pcmp=1; % Component for plotting of u=[u,v]
huclean(p); % Set up figures in nice way

p.sw.sfem=-1; % Use OOPDE settings. Otherwise 0/1 for full/preassembled FEM settings
p.nc.neq=2; % Number of equations

p.fuha.sG=@sG; % Functions to be solved
p.sw.jac=0; % Use numeric Jacobian.

% Space setup and number of elements to be used, space is [-lx,lx]X[-ly,ly].
p.dim=2; % Spatial dimensions
A=2.5;

sw.sym=2;
pde=EllipsePDE(A,B,100,0.3); % PDE setup
p.pdeo=pde;p.np=pde.grid.nPoints;p.nu=2*p.np;p.nt=pde.grid.nElements;
p.sw.verb=0; % Verbosity of calculation
p.sol.xi=1/p.nu; % Normalisation weights

p.sw.bifcheck=2; % Calculate bifurcations using eigenvalues, 0 turn off, 1 use LU decomposition
p.sw.spcalc=1; % Calculate eigenvalues (0/1 Eigenvalue computations off/on)
p.nc.neig=min(50,p.np); % Number of eigenvalues to calculate
% p.nc.neig=min(20,p.np); % Number of eigenvalues to calculate
p.nc.tol=1e-11; % Tolerance for finding branch
p.nc.nsteps=1000; % Number of continuation steps

p.sol.ds=1e-4; % Arclength step
p.nc.dsmin=1e-6; % Arclength minimum steplength
p.nc.dsmax=0.01; % Arclength max steplength
p.nc.dlammax=0.1; % Maximum step in the parameter
p.nc.lammax=1.08; % Max parameter value
p.nc.lammin=1; % Min parameter value
p.nc.mu2=0.01; % Threshold for bifcheck=2
p.file.smod=10; % Output every n steps

p.sw.para=1; % Continuation variable 1: automatic switching via ? <> p.nc.lamdtol (0: natural parameter.; 2: arclength).
p.plot.pstyle=2; %0 only plot the FEM mesh; 1 mesh-plot; 2 density-plot

u=(par(4)+par(5))*ones(p.np,1); v=par(5)/(par(4)+par(5))^2*ones(p.np,1);  %ICs
p.u=[u; v; par']; % Complete description

%NEED TO CHANGE oosetfemops
p=oosetfemops(p);
p=setfn(p,Fname,0.01);

% p=findbif(p,1);
p=cont(p);
%%
p=swibra(Fname,'bpt1',[Fname,'bp1'],0.01);
p=cont(p,1000);




function r=sG(p,u)
n=p.np;


par=u(p.nu+1:end);
L=par(1);
Du=par(2);Dv=par(3);
alpha=par(4);
beta=par(5);


U=u(1:n); V=u(n+1:2*n);

f1=alpha-U+U.^2.*V;
f2=beta-U.^2.*V;
f=[f1;f2]; % semilin.nonlinearity


Kx=p.mat.Kx;
Ky=p.mat.Ky;
r=[Du*(Kx/L^2+Ky) 0*Kx; 0*Kx Dv*(Kx/L^2+Ky)]*[U;V]-p.mat.M*f; % putting rhs together
end

function out=outfn(p,u)
% output to bifurcation diagram function
% u=u(1:p.nu);n=p.np;u1=u(1:n); u2=u(n+1:2*n); % Seperate parameters and variables
out=[u(p.nu+1:end); % parameters
    max(u(1:p.np));
    min(u(1:p.np))];
end
