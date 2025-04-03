%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a simple example of angular Mathieu function Spm computation
% at given:  order n, nmax, and different values of q
% [( nmax >= n )& (nmax <= 25)]
%
% Results are shown in Fig.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    E. Cojocaru, revised November 2008
%    Observations, suggestions, and recommendations are welcome at e-mail:
%    ecojocaru@theory.nipne.ro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters 
n=1;                   % single value of order n
vq=1:5;                % 5 values of elliptical parameter q

nmax=1;                % nmax ( >= n and <= 25)  
vv=0:pi/20:2*pi;       % 41 values of angle v in radians 

% Specify the function code k:
for k=1:4
  
    my=[];
    
    for kq=1:5
        
        q=vq(kq);             % take a value of q 
        
[va,mv,vt]=eig_Spm(k,q);      % compute characteristic values and
                              % expansion coefficients 

y=[];

for kv=1:length(vv);
v=vv(kv);                     % take a value of angle v
vy=Spm(k,v,mv,nmax);          % compute Spm for nmax orders; size(vy)=[1 1]

yn=vy(n);                     % extract values of Spm corresponding to 
                              % order n; size(yn)=[1 1];

y=[y; yn];                    % column vector of Spm values at different
                              % angles v, at order n, and a value q;
                              % size(y)=[41 1]
                              
end


my=[my y];                    % matrix of Spm values at different
                              % values of v and q; size(my)=[41 5]  

end

% plot results

subplot(2,2,k);
plot(vv/pi,my(:,1),'k.-',vv/pi,my(:,2),'r.-',vv/pi,my(:,3),'b.-', ...
     vv/pi,my(:,4),'c.-',vv/pi,my(:,5),'m.-'); hold on
 
xlabel('v/pi')
ylabel('Spm')
title([' n=1; q=1:5; KF=',num2str(k)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
