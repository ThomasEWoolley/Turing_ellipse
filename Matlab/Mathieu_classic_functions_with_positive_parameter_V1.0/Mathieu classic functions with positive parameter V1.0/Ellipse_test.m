ccc
% Input parameters
q=1;
vv=linspace(0,2*pi);
uu=linspace(0,3);
[uvec,vvec]=meshgrid(uu,vv);
x=cosh(uvec).*cos(vvec);
y=sinh(uvec).*sin(vvec);
i=1;
for k=[1 2]
    for n=[1 2]

    if i==1
        uu=linspace(0,1.69);
    elseif i==2
        uu=linspace(0,1.09);
    elseif i==3
        uu=linspace(0,0.8);
    elseif i==4
        uu=linspace(0,1.15);
    end

    % if i==1
    %     uu=linspace(0,1.33);
    % elseif i==2
    %     uu=linspace(0,1.6);
    % elseif i==3
    %     uu=linspace(0,1.21);
    % elseif i==4
    %     uu=linspace(0,0.81);
    % end


vv=linspace(0,2*pi);
[uvec,vvec]=meshgrid(uu,vv);
x=cosh(uvec).*cos(vvec);
y=sinh(uvec).*sin(vvec);


        my=mathieu(k,n,q,vv,'what','what');
        ny=mmathieu1(k,n,q,uu);

        [R,S]=meshgrid(ny,my);

        % plot results
        figure(1)
        subplot(2,2,i);
        
        plot(vv/pi,my)
        % set(gca,'XLim',[0 0.5],'YLim',[-1 1])
        xlabel('v/pi')
        ylabel('Spm')

        figure(2)
        subplot(2,2,i);
        plot(uu,ny)
        % set(gca,'XLim',[0 0.5],'YLim',[-1 1])
        xlabel('u')
        ylabel('Rpm')

        figure(3)

        subplot(2,2,i);
        pcolor(x,y,R.*S)
        shading interp
        % set(gca,'XLim',[0 0.5],'YLim',[-1 1])
        xlabel('x')
        ylabel('y')
        i=i+1;
    end
end
%
%
% % Specify the function code k:
% for k=1:4
%
%     if k==1
%         uu=linspace(0,1.35);
%     elseif k==2
%         uu=linspace(0,1.65);
%     elseif k==3
%         uu=linspace(0,1.12);
%     elseif k==4
%         uu=linspace(0,0.76);
%     end
%
%
%     [uvec,vvec]=meshgrid(uu,vv);
%     x=cosh(uvec).*cos(vvec);
%     y=sinh(uvec).*sin(vvec);
%
%
%     [~,mv,~]=eig_Spm(k,q);      % computes characteristic values and
%     % expansion coefficients
%     my=[];
%     ny=[];
%
%     for i=1:length(vv);
%
%         vy=Spm(k,vv(i),mv,n);
%         uy=Jpm(k,uu(i),q,mv,n);
%
%         my=[my; vy(n)];                  % matrix of Spm values at different values of
%         ny=[ny; uy(n)];                  % matrix of Spm values at different values of
%         % v and n;  size(my)=[51 3]
%
%     end
%
%     [R,S]=meshgrid(ny,my);
%
%     % plot results
%     figure(1)
%     subplot(2,2,k);
%     plot(vv/pi,my)
%     % set(gca,'XLim',[0 0.5],'YLim',[-1 1])
%     xlabel('v/pi')
%     ylabel('Spm')
%
%     figure(2)
%     subplot(2,2,k);
%     plot(uu,ny)
%     % set(gca,'XLim',[0 0.5],'YLim',[-1 1])
%     xlabel('u')
%     ylabel('Rpm')
%
%     figure(3)
%
%     subplot(2,2,k);
%     pcolor(x,y,R.*S)
%     shading interp
%     % set(gca,'XLim',[0 0.5],'YLim',[-1 1])
%     xlabel('x')
%     ylabel('y')
%
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%