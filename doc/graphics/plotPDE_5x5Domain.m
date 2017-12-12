clear
close all

N = 5;
x = 1:N;
y = 1:N;

%%% Draw base grid
fig = figure;
hold on;
[XX, YY] = meshgrid(x,y);
mPlot = mesh(XX,YY,zeros(size(XX)));
mPlot.EdgeColor = [0 0 0];
mPlot.Marker = 'o';
mPlot.MarkerSize = 12;
mPlot.MarkerFaceColor = [0 0 0];
mPlot.MarkerFaceColor = [1 1 1];

%%% Place Grid Indices
kk = 0;
for ii = 1:N
    for jj = 1:N
        kk = kk+1;
        IND{kk} = text(XX(kk),YY(kk),num2str(kk));
        IND{kk}.HorizontalAlignment = 'center';
        IND{kk}.FontName = 'FixedWidth';
        IND{kk}.FontWeight = 'Bold';
    end
    
end

%%% Place index counters around the outside
% X-Dir counters
idx = 0;
for ii = 1:N
    % NORTH
    idx = idx+1;
    str = "i = " + num2str(ii);
    IDX{idx} = text(XX(end,ii),YY(end,ii)+0.25,str);
    IDX{idx}.HorizontalAlignment = 'center';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % SOUTH
    idx = idx+1;
    str = "i = " + num2str(ii);
    IDX{idx} = text(XX(1,ii),YY(1,ii)-0.25,str);
    IDX{idx}.HorizontalAlignment = 'center';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
end

% Y-Dir counters
idx = 0;
for jj = 1:N
    % WEST
    idx = idx+1;
    str = "j = " + num2str(jj);
    IDX{idx} = text(XX(jj,1)-0.25,YY(jj,1),str);
    IDX{idx}.HorizontalAlignment = 'right';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % EAST
    idx = idx+1;
    str = "j = " + num2str(jj);
    IDX{idx} = text(XX(jj,end)+0.25,YY(jj,end),str);
    IDX{idx}.HorizontalAlignment = 'left';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
end

axis ij
axis equal
axis off

saveas(fig,'PDE_5x5Domain','svg')