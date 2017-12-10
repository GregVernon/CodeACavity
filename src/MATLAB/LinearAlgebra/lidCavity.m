function lidCavity(N, tmax, method, tol, max_tstep, plotFlow)

NY = N;
NX = N;

xmin = 0;
xmax = 1;
x = linspace(xmin,xmax,NX);
ymin = 0;
ymax = 1;
y = linspace(ymin,ymax,NY);
[xx,yy] = meshgrid(x,y);

dx=x(2)-x(1);
dy=y(2)-y(1);

Visc=.01;

dt_max = tmax / 10000;

PSI=zeros(NY,NX);
OMEGA=zeros(NY,NX);
OMEGA([NY NY-1],2:NX-1)= -2.0/dy;

U = zeros(NY,NX);
V = zeros(NY,NX);
VELOCITY = zeros(NY,NX);
cREYNOLDS = zeros(NY,NX);

FDM = assembleCoeffMatrix(dx, dy, NX,NY);

fixedPointMethods = {'Jacobi','Weighted Jacobi','Mapped Jacobi','Richardson','Mapped Richardson','AlternatingAndersonRichardson'};
if strcmpi(method,'Decomposition')
    FDM = decomposition(FDM);
elseif any(ismember(fixedPointMethods,method))
    A = FDM;
    clearvars FDM
    FDM.A = A;
    dVal = spdiags(A,0);
    I = speye(size(A));
    [iRow,iCol] = find(I);
    D = sparse(iCol,iRow,dVal,size(A,1),size(A,2));
    FDM.D = D;
    FDM.R = triu(A,1) + tril(A,-1);
    FDM.iD = inv(FDM.D);
end

initPlot = true;
tic;
t=0.0;
pIter = 0;
tstep = 0;
while t < tmax && tstep < max_tstep% start the time integration
    pIter = pIter+1;
    tstep = tstep + 1;
    % Compute timestep
    dt = delta_t(VELOCITY, cREYNOLDS, dt_max, dx,dy);
    % Compute streamfunction
    PSI = computePSI(PSI, OMEGA, FDM, method,tol, dx, dy, NX, NY);
    % Apply vorticity boundary conditions
    OMEGA = applyBC_OMEGA(PSI, OMEGA, dx, dy, NX, NY);
    % Compute Velocity
    [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, dx, dy, NX, NY);
    % Compute vorticity
    OMEGA = computeOMEGA(PSI, OMEGA, U, V, dt, dx, dy, NX, NY, Visc);
    % Increment time value by timestep
    t=t+dt;
    %% plot
    if pIter == 1e2
        pIter = 0;
        disp(['Time: ' num2str(t)])
        if plotFlow == true
            if initPlot == true
                initPlot = false;
                figure;
                subplot(131);
                [~,pltOMEGA]=contourf(xx,yy,OMEGA); 
                axis('square'); 
                colormap(jet);
                caxis([-max(abs(OMEGA(:))) max(abs(OMEGA(:)))]); 
                colorbar('southoutside')% plot vorticity
                
                subplot(132);
                [~,pltPSI]=contourf(xx,yy,PSI); 
                axis('square'); 
                colormap(jet);
                caxis([-max(abs(PSI(:))) max(abs(PSI(:)))]);
                colorbar('southoutside')% streamfunction
                
                subplot(133);
                hold on
                pltVEL=quiver(xx,yy,U,V); 
                [~,conVEL]=contour(xx,yy,VELOCITY);
                axis('square'); 
                axis([xmin xmax ymin ymax]);
                colorbar('southoutside')
                drawnow % streamfunction
            else
                pltOMEGA.ZData=OMEGA;
                pltPSI.ZData=PSI;
                pltVEL.UData=U; pltVEL.VData=V; conVEL.ZData=VELOCITY; drawnow
            end
        end
    end
end