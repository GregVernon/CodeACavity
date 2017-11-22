%%% Draw domain
clear
close all
fig = figure; hold on
x = [1 0 0    0 1 1      1  0];
y = [1 1 0.985 0 0 0.985 1  1];
dPlot = plot(x,y); drawnow
dPlot.LineWidth = 4;
dPlot.AlignVertexCenters = 'on';
dPlot.LineJoin = 'miter';
cd = uint8([255 255 0 0 0  0 255 255; ...
             0    0 0 0 0  0   0   0 ; ...
             0    0 0 0 0  0   0   0 ;...
             1    1 1 1 1  1   1   1]);
set(dPlot.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
drawnow

% Top Wall
t{1} = text(0.5,0.9,["$$\frac{\partial\Psi}{\partial y} = U$$ ; $$\omega = -\frac{\partial^{2}\Psi}{\partial y^{2}}$$"],'interpreter','latex');
t{1}.HorizontalAlignment = 'center';
t{1}.FontSize = 16;
t{1}.Color = [1 0 0];
% Bottom Wall
t{2} = text(0.5,0.1,'$$\frac{\partial\Psi}{\partial y} = 0$$ ; $$\omega = -\frac{\partial^{2}\Psi}{\partial y^{2}}$$','interpreter','latex');
t{2}.HorizontalAlignment = 'center';
t{2}.FontSize = 16;
% Left Wall
t{3} = text(0.1,0.5,'$$\frac{\partial\Psi}{\partial x} = 0$$  ; $$\omega = -\frac{\partial^{2}\Psi}{\partial x^{2}}$$','interpreter','latex');
t{3}.HorizontalAlignment = 'center';
t{3}.Rotation = 90;
t{3}.FontSize = 16;
% Right Wall
t{4} = text(0.9,0.5,'$$\frac{\partial\Psi}{\partial x} = 0$$  ; $$\omega = -\frac{\partial^{2}\Psi}{\partial x^{2}}$$','interpreter','latex');
t{4}.HorizontalAlignment = 'center';
t{4}.Rotation = 90;
t{4}.FontSize = 16;

axis equal
axis off
saveas(fig,'Domain_Definition','svg')