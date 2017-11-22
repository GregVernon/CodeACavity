clear
close all

N = 3;
x{1} = 1:N;
x{2} = 1:N;
x{3} = N+1 : N+N;
x{4} = N+1 : N+N;

y{1} = 1:N;
y{2} = N+1 : N+N;
y{3} = 1:N;
y{4} = N+1 : N+N;

%%% Draw base grid
fig = figure;
hold on;
for ii = 1:4
    [XX{ii}, YY{ii}] = meshgrid(x{ii},y{ii});
    mPlot{ii} = mesh(XX{ii},YY{ii},zeros(size(XX{ii})));
    mPlot{ii}.EdgeColor = [0 0 0];
    mPlot{ii}.Marker = 'o';
    mPlot{ii}.MarkerSize = 12;
    mPlot{ii}.MarkerFaceColor = [0 0 0];
    mPlot{ii}.MarkerFaceColor = [1 1 1];
end

%%% Draw connectors
% Connect SW to SE
jj = 0;
for ii = 1:N
    jj = jj+1;
    cPlot{jj} = plot([XX{1}(ii,end) XX{3}(ii,1)], [YY{1}(ii,end) YY{3}(ii,1)]);
    cPlot{jj}.LineStyle = ':';
    cPlot{jj}.Color = [0 0 0];
    uistack(cPlot{jj},'down',1e6)
end

% Connect NW to NE
for ii = 1:N
    jj = jj+1;
    cPlot{jj} = plot([XX{2}(ii,end) XX{4}(ii,1)], [YY{2}(ii,end) YY{4}(ii,1)]);
    cPlot{jj}.LineStyle = ':';
    cPlot{jj}.Color = [0 0 0];
    uistack(cPlot{jj},'down',1e6)
end

% Connect SW to NW
for ii = 1:N
    jj = jj+1;
    cPlot{jj} = plot([XX{1}(end,ii) XX{2}(1,ii)], [YY{1}(end, ii) YY{2}(1,ii)]);
    cPlot{jj}.LineStyle = ':';
    cPlot{jj}.Color = [0 0 0];
    uistack(cPlot{jj},'down',1e6)
end

% Connect SE to NE
for ii = 1:N
    jj = jj+1;
    cPlot{jj} = plot([XX{3}(end,ii) XX{4}(1,ii)], [YY{3}(end,ii) YY{4}(1,ii)]);
    cPlot{jj}.LineStyle = ':';
    cPlot{jj}.Color = [0 0 0];
    uistack(cPlot{jj},'down',1e6)
end

%%% Place Alphabetical Letters in the grid
% Start Top left, move left-to-right first, then top-to-bottom
upABC = uint8('A'):uint8('Z');
loABC = uint8('a'):uint8('z');
abc = upABC(1)-1;

row = N;
for jj = 1:N-1
    row = row - 1;
    % NW Quadrant
    col = 1;
    for ii = 1:N-1
        abc = abc+1;
        if abc > max(upABC) && abc < min(loABC)
            abc = loABC(1);
        end
        col = col + 1;
        
        ABC{abc} = text(XX{2}(row,col),YY{2}(row,col),char(abc));
        ABC{abc}.HorizontalAlignment = 'center';
        ABC{abc}.FontName = 'FixedWidth';
        ABC{abc}.FontWeight = 'Bold';
    end
    
    % NE Quadrant
    col = 0;
    for ii = 1:N-1
        abc = abc+1;
        if abc > max(upABC) && abc < min(loABC)
            abc = loABC(1);
        end
        col = col + 1;
        
        ABC{abc} = text(XX{4}(row,col),YY{4}(row,col),char(abc));
        ABC{abc}.HorizontalAlignment = 'center';
        ABC{abc}.FontName = 'FixedWidth';
        ABC{abc}.FontWeight = 'Bold';
    end
end

row = N+1;
for jj = 1:N-1
    row = row - 1;
    % SW Quadrant
    col = 1;
    for ii = 1:N-1
        abc = abc+1;
        if abc > max(upABC) && abc < min(loABC)
            abc = loABC(1);
        end
        col = col + 1;
        
        ABC{abc} = text(XX{1}(row,col),YY{1}(row,col),char(abc));
        ABC{abc}.HorizontalAlignment = 'center';
        ABC{abc}.FontName = 'FixedWidth';
        ABC{abc}.FontWeight = 'Bold';
    end
    
    % SE Quadrant
    col = 0;
    for ii = 1:N-1
        abc = abc+1;
        if abc > max(upABC) && abc < min(loABC)
            abc = loABC(1);
        end
        col = col + 1;
        
        ABC{abc} = text(XX{3}(row,col),YY{3}(row,col),char(abc));
        ABC{abc}.HorizontalAlignment = 'center';
        ABC{abc}.FontName = 'FixedWidth';
        ABC{abc}.FontWeight = 'Bold';
    end
end

%%% Place index counters around the outside
% X-Dir counters
idx = 0;
for ii = 1:N
    % NW Quadrant
    idx = idx+1;
    str = "i = " + num2str(ii);
    IDX{idx} = text(XX{2}(end,ii),YY{2}(end,ii)+0.25,str);
    IDX{idx}.HorizontalAlignment = 'center';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % SW Quadrant
    idx = idx+1;
    str = "i = " + num2str(ii);
    IDX{idx} = text(XX{1}(1,ii),YY{1}(1,ii)-0.25,str);
    IDX{idx}.HorizontalAlignment = 'center';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % NE Quadrant
    idx = idx+1;
    if ii == N
        str = "i = NX";
    else
        str = "i = NX-" + num2str(N-ii);
    end
    IDX{idx} = text(XX{4}(end,ii),YY{4}(end,ii)+0.25,str);
    IDX{idx}.HorizontalAlignment = 'center';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % SE Quadrant
    idx = idx+1;
    if ii == N
        str = "i = NX";
    else
        str = "i = NX-" + num2str(N-ii);
    end
    IDX{idx} = text(XX{3}(1,ii),YY{3}(1,ii)-0.25,str);
    IDX{idx}.HorizontalAlignment = 'center';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
end

% Y-Dir counters
idx = 0;
for jj = 1:N
    % NW Quadrant
    idx = idx+1;
    str = "j = " + num2str((N+1)-jj);
    IDX{idx} = text(XX{2}(jj,1)-0.25,YY{2}(jj,1),str);
    IDX{idx}.HorizontalAlignment = 'right';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % SW Quadrant
    idx = idx+1;
    if jj == 1
        str = "j = NY";
    else
        str = "j = NY-" + num2str(jj-1);
    end
    IDX{idx} = text(XX{1}(jj,1)-0.25,YY{1}(jj,1),str);
    IDX{idx}.HorizontalAlignment = 'right';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % NE Quadrant
    idx = idx+1;
    str = "j = " + num2str((N+1)-jj);
    IDX{idx} = text(XX{4}(jj,end)+0.25,YY{4}(jj,end),str);
    IDX{idx}.HorizontalAlignment = 'left';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
    
    % SE Quadrant
    idx = idx+1;
    if jj == 1
        str = "j = NY";
    else
        str = "j = NY-" + num2str(jj-1);
    end
    IDX{idx} = text(XX{3}(jj,end)+0.25,YY{3}(jj,end),str);
    IDX{idx}.HorizontalAlignment = 'left';
    IDX{idx}.FontName = 'FixedWidth';
    IDX{idx}.FontWeight = 'Bold';
    IDX{idx}.FontSize = 9;
end

axis equal
axis off

saveas(fig,'PDE_IntDomain','svg')