%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a simple example of angular Mathieu function Spm computation
% at given:  q > 0, nmax, and different values of order n
% [( nmax >= max(n) ) & (nmax <= 25)]
%
% Results are shown in Fig.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters 
q=1;                      % elliptical parameter 
n=3:5;                    % 3 different values of order n
nmax=5;                   % nmax ( >= max(n) and <= 25)
vv=0:pi/100:pi/2;         % 51 values of angle v in radians

% Specify the function code k:
for k=1:4
    
[va,mv,vt]=eig_Spm(k,q);      % computes characteristic values and
                              % expansion coefficients 
                              
my=[];                    

for kv=1:length(vv);
    
 v=vv(kv);                    % take a value of angle v
 
vy=Spm(k,v,mv,nmax);          % compute Spm for nmax orders; size(vy)=[1 5]
yn=vy(n);                     % extract values of Spm corresponding to 
                              % orders n=3:5; size(yn)=[1 3]
                              
my=[my; yn];                  % matrix of Spm values at different values of
                              % v and n;  size(my)=[51 3]                         
                          
end

% plot results

subplot(2,2,k);

plot(vv/pi,my(:,1),'k.-',vv/pi,my(:,2),'r.-',vv/pi,my(:,3),'b.-'); hold on

set(gca,'XLim',[0 0.5],'YLim',[-1 1]) 
xlabel('v/pi')
ylabel('Spm')
title([' q=1; n=3:5; KF=',num2str(k)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%