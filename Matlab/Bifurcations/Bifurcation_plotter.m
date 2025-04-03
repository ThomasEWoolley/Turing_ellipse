clear all
close all
clc
figure('Units','normalized','Position',[1 0 1 1/3])
subplot(1,5,1:2)
hold on
cols=[1 0 0;
    0 0.5 0;
    0 0 1;
    200/256, 200/256, 0
    203/256, 195/256, 227/256
    0 0 0];


nos=[0.5, 1, 1.5,2];
l=1;
for j=nos
    Name=['./Ellipsehom_',num2str(j)];
    maxNumber=NoFinder(Name);

    load(['./Ellipsehom_',num2str(j),'/pt',num2str(maxNumber),'.mat'])
    n1=1;
    n2=length(p.branch(8,:));
    % pb=plot(2.5*p.branch(4,n1:n2),p.branch(end-1,n1:n2),'linewidth',3); %Plots branch max
    % plot(2.5*p.branch(4,n1:n2),p.branch(end,n1:n2),'linewidth',3) %Plots branch min

    for i=n1:n2-1 % loop from first-point to last point
        lw=1;
        if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
        if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end

        plot(2.5*p.branch(4,i:i+1),p.branch(end-1,i:i+1),'Linewidth',lw,'Color',cols(l,:));
        plot(2.5*p.branch(4,i:i+1),p.branch(end,i:i+1),'Linewidth',lw,'Color',cols(l,:));

    end

    Name=['./Ellipsehom_',num2str(j),'bp1'];
    maxNumber=NoFinder(Name);

    load(['./Ellipsehom_',num2str(j),'bp1/pt',num2str(maxNumber),'.mat'])

    n1=1;
    n2=length(p.branch(8,:));

    pb=plot(2.5*p.branch(4,n1:n2),p.branch(end-1,n1:n2),'linewidth',3,'Color',cols(l,:)); %Plots branch max
    plot(2.5*p.branch(4,n1:n2),p.branch(end,n1:n2),'linewidth',3,'Color',cols(l,:)) %Plots branch min
    Bifpoint(l)=2.5*p.branch(4,1);
    l=l+1;
end
j=2.6038;
Name=['./Ellipsehom_circle'];
maxNumber=NoFinder(Name);

load(['./Ellipsehom_circle/pt',num2str(maxNumber),'.mat'])
n1=1;
n2=length(p.branch(8,:));
% pb=plot(2.5*p.branch(4,n1:n2),p.branch(end-1,n1:n2),'linewidth',3); %Plots branch max
% plot(2.5*p.branch(4,n1:n2),p.branch(end,n1:n2),'linewidth',3) %Plots branch min

for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end

    plot(p.branch(4,i:i+1),p.branch(end-1,i:i+1),'Linewidth',lw,'Color',cols(l,:));
    plot(p.branch(4,i:i+1),p.branch(end,i:i+1),'Linewidth',lw,'Color',cols(l,:));

end

Name=['./Ellipsehom_circlebp1'];
maxNumber=NoFinder(Name);

load(['./Ellipsehom_circlebp1/pt',num2str(maxNumber),'.mat'])

n1=1;
n2=length(p.branch(8,:));

pb=plot(p.branch(4,n1:n2),p.branch(end-1,n1:n2),'linewidth',3,'Color',cols(l,:)); %Plots branch max
plot(p.branch(4,n1:n2),p.branch(end,n1:n2),'linewidth',3,'Color',cols(l,:)) %Plots branch min
Bifpoint(l)=p.branch(4,1);
l=l+1;



Name=['./Ellipsehom_line'];
maxNumber=NoFinder(Name);

load(['./Ellipsehom_line/pt',num2str(maxNumber),'.mat'])
n1=1;
n2=length(p.branch(8,:));
% pb=plot(2.5*p.branch(4,n1:n2),p.branch(end-1,n1:n2),'linewidth',3); %Plots branch max
% plot(2.5*p.branch(4,n1:n2),p.branch(end,n1:n2),'linewidth',3) %Plots branch min

for i=n1:n2-1 % loop from first-point to last point
    lw=1;
    if(~ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i)  ==0); lw=3;  end
    if( ismember(p.branch(2,i),[-2 -1 1 2]) && p.branch(3,i+1)==0); lw=3;  end

    plot(p.branch(4,i:i+1),p.branch(end-1,i:i+1),'Linewidth',lw,'Color',cols(l,:));
    plot(p.branch(4,i:i+1),p.branch(end,i:i+1),'Linewidth',lw,'Color',cols(l,:));

end

Name=['./Ellipsehom_linebp1'];
maxNumber=NoFinder(Name);
load(['./Ellipsehom_linebp1/pt',num2str(maxNumber),'.mat'])

n1=1;
n2=length(p.branch(8,:));

pb=plot(p.branch(4,n1:n2),p.branch(end-1,n1:n2),'linewidth',3,'Color',cols(l,:)); %Plots branch max
plot(p.branch(4,n1:n2),p.branch(end,n1:n2),'linewidth',3,'Color',cols(l,:)) %Plots branch min
Bifpoint(l)=p.branch(4,1);
l=l+1;



% plot(Bifpoint,1,'s')
for i=1:4
    cols(i,:)
    plot(Bifpoint(i),1,'d','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',10)
end
i=5;
plot(Bifpoint(i),1,'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',10)
i=6;
plot(Bifpoint(i),1,'s','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',10)

p = zeros(1,6);
for i=1:6
    p(i)=plot([0 1],[0.1 0.1],'color',cols(i,:),'LineWidth',3)
end

legend(p,'$B=0.5$','$B=1$','$B=1.5$','$B=2$','Circle','Line')
xlim([2.2 2.7])
xlabel('$A$')
ylabel('Steady state amplitudes, $u$')
set(gca,'fontsize',12)
export_fig(['../../Pictures/Bifs.png'],'-r300')
% set(h,'FontSize',15);
%%
figure('Units','normalized','Position',[1 0 1 1/3])
subplot(1,5,3)
p=loadp('Ellipsehom_0.5bp1');
pdeplot(p.branch(4,end)*p.pdeo.grid.p,p.pdeo.grid.e,p.pdeo.grid.t,'xydata',p.u(1:p.np));
axis equal
colorbar off
colorbar
colormap parula
xlabel('$x$')
ylabel('$y$')
set(gca,'fontsize',12)
export_fig(['../../Pictures/Ellipse.png'],'-r300')

figure('Units','normalized','Position',[1 0 1 1/3])
subplot(1,5,4)
p=loadp('Ellipsehom_circlebp1');
pdeplot(p.branch(4,end)*p.pdeo.grid.p,p.pdeo.grid.e,p.pdeo.grid.t,'xydata',p.u(1:p.np));
axis equal
colorbar off
colorbar
xlabel('$x$')
ylabel('$y$')
colormap parula
set(gca,'fontsize',12)
export_fig(['../../Pictures/Circle.png'],'-r300')

figure('Units','normalized','Position',[1 0 1 1/3])
subplot(1,5,5) %Change oosetfemops
p=loadp('Ellipsehom_linebp1');
% plot(p.branch(4,end)*p.pdeo.grid.p,p.u(1:p.np),'k','LineWidth',3);

t = p.branch(4,end)*p.pdeo.grid.p; % Time variable
y = p.u(1:p.np); % Function

hold on;
cmap = parula(256); 
y_min = min(y);
y_max = max(y);

% Map y values to the colormap indices without normalisation
color_indices = round(((y - y_min) / (y_max - y_min)) * (size(cmap, 1) - 1)) + 1;

% Plot each segment with the corresponding colour
for i = 1:length(t)-1
    plot(t(i:i+1), y(i:i+1), 'Color', cmap(color_indices(i), :), 'LineWidth', 3);
end

% Add colorbar to indicate mapping
colorbar
caxis([y_min y_max]); % Set colorbar limits to match y range
axis tight
xlabel('$x$')
ylabel('$u$')
set(gca,'fontsize',12)





export_fig(['../../Pictures/Line.png'],'-r300')

function maxNumber=NoFinder(folderPath)

% Get a list of all .mat files in the folder
fileList = dir(fullfile(folderPath, 'pt*.mat'));

% Initialise an array to store the extracted numbers
numbers = [];

% Loop through the file names and extract numbers
for i = 1:length(fileList)
    fileName = fileList(i).name;
    match = regexp(fileName, 'pt(\d+)\.mat$', 'tokens');
    if ~isempty(match)
        numbers(end+1) = str2double(match{1}{1}); % Convert to number and store
    end
end


maxNumber = max(numbers);

end