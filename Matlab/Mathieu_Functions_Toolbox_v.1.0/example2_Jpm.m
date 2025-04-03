%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a simple example of radial Mathieu function Jpm computation
% at given:  n, nmax, and different values of elliptical parameter q
% [( nmax >= max(n) ) & (nmax <= 25)]
%
% Results are shown in Fig.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters 
n=1;                       % single value of n
vq=1:5;                    % 5 values of elliptical parameter q

nmax=1;                    % nmax ( >= max(n) and  <= 25)
vu=0:5e-02:2.5;            % 51 values of coordinate u 

% Specify the function code k:
for k=1:4
  
    my=[];
    
    for kq=1:5
        q=vq(kq);           % take a value of q
[va,mv,vt]=eig_Spm(k,q);    % compute characteristic values and
                            % expansion coefficients

y=[];

for ku=1:length(vu);
u=vu(ku);                   % take a value of u
vy=Jpm(k,u,q,mv,nmax);      % compute Jpm for nmax orders; size(vy)=[1 1]

yn=vy(n);                   % extract values of Jpm corresponding to n;
                            % size(yn)=[1 1];

y=[y; yn];                  % column vector of Jpm values at different
                            % values of u, at order n, and a value q;
                            % size(y)=[51 1] 

end

my=[my y];                  % matrix of Jpm values at different
                            % values of u and q; size(my)=[51 5]  

end

% plot results

subplot(2,2,k);
plot(vu,my(:,1),'k.-',vu,my(:,2),'r.-',vu,my(:,3),'b.-', ...
     vu,my(:,4),'c.-',vu,my(:,5),'m.-'); hold on
 
set(gca,'XLim',[0 2.5])  
xlabel('u')
ylabel('Jpm')
title([' n=1; q=1:5; KF=',num2str(k)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
