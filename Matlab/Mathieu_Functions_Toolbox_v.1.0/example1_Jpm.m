%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a simple example of radial Mathieu function Jpm computation
% at given:  q > 0, nmax, and different values of order n
% [( nmax >= max(n) ) & (nmax <= 25)]
%
% Results are shown in Fig.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters 
q=1;                       % single value of elliptical parameter q
n=[3 6 9];                 % 3 values of order n 
nmax=9;                    % nmax ( >= max(n) and  <= 25)
vu=0:5e-02:2.5;            % 51 values of coordinate u 

% Specify the function code k:
for k=1:4
    
[va,mv,vt]=eig_Spm(k,q);      % compute characteristic values and
                              % expansion coefficients

my=[];

for ku=1:length(vu);
 u=vu(ku);                    % take a value of u
 
vy=Jpm(k,u,q,mv,nmax);        % compute Jpm for nmax orders; size(vy)=[1 9]

yn=vy(n);                     % extract values of Jpm corresponding to 
                              % orders n=3,6,9; size(yn)=[1 3]    

my=[my; yn];                  % matrix of Jpm values at different values of
                              % u and n; size(my)=[51 3]

end

% plot results

subplot(2,2,k);

plot(vu,my(:,1),'k.-',vu,my(:,2),'r.-',vu,my(:,3),'b.-'); hold on

set(gca,'XLim',[0 2.5]) 
xlabel('u')
ylabel('Jpm')
title([' q=1; n=3,6,9; KF=',num2str(k)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
