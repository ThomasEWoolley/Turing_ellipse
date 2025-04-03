function p=oosetfemops(p)
gr=p.pdeo.grid; fem=p.pdeo.fem;n=p.np;

% [K,M,~]=fem.assema(gr,1,1,1); %Isometric diffusion

% L=p.u(p.nu+1);
% c1=1/L^2; c2=0; c3=0; c4=1; % diffusion coeffients
% cc=[[c1 c2]; [c3 c4]]; % stored as a "2 x (2*np) matrix", see grid2D.aCoefficientsMpt
% [K,M,~]=fem.assema(gr,cc,1,1); %Anisometric diffusion
[K,M,~]=fem.assema(gr,1,1,1);
p.mat.M=[M 0*M; 0*M M];
[Kx,~,~]=fem.assema(gr,[ones(n,1);0*ones(n,1)],1,1);
[Ky,~,~]=fem.assema(gr,[0*ones(n,1); ones(n,1)],1,1); 
p.mat.Kx=Kx;
p.mat.Ky=Ky;
E=center2PointMatrix(gr); % to map element differentiation matrices to nodal ones
p.mat.p2c=point2CenterMatrix(gr); % to interpolate from nodes to element centers
end
