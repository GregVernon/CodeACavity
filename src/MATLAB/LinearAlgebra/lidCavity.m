function lidCavity(N, tmax, method, tol, max_tstep, plotFlow)

NY = N;
NX = N;
h=1.0/(NY-1);

Visc=.01;

dt_max = tmax / 10000;

PSI=zeros(NY,NX);
OMEGA=zeros(NY,NX);
jj = 2:NY-1;
OMEGA(jj,[NX-1 NX])= -2.0/h;

VELOCITY = zeros(NY,NX);
cREYNOLDS = zeros(NY,NX);

FDM = assembleCoeffMatrix(NX,NY);

if strcmpi(method,'Decomposition')
    FDM = decomposition(FDM);
elseif strcmpi(method,'Jacobi')
    A = FDM;
    clearvars FDM
    FDM.A = A;
    FDM.D = speye(size(A)) .* spdiags(A,0);
    FDM.R = triu(A,1) + tril(A,-1);
    FDM.iD = inv(FDM.D);
end

tic;
t=0.0;
pIter = 0;
tstep = 0;
while t < tmax && tstep < max_tstep% start the time integration
    pIter = pIter+1;
    tstep = tstep + 1;
    % Compute timestep
    dt = delta_t(VELOCITY, cREYNOLDS, dt_max, h);
    % Compute streamfunction
    PSI = computePSI(PSI, OMEGA, FDM, method,tol, h, NY, NX);
    % Apply vorticity boundary conditions
    OMEGA = applyBC_OMEGA(PSI, OMEGA, h, NY, NX);
    % Compute vorticity
    OMEGA = computeOMEGA(PSI, OMEGA, dt, h, NY, NX, Visc);
    % Compute Velocity
    [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, h, NY, NX);
    % Increment time value by timestep
    t=t+dt;
    %% plot
    if pIter == 1e3
        pIter = 0;
        disp(['Time: ' num2str(t)])
        if plotFlow == true
            subplot(131), contourf(rot90(fliplr(OMEGA))), axis('square'); colorbar% plot vorticity
            subplot(132), contourf(rot90(fliplr(PSI))), axis('square'); colorbar% streamfunction
            subplot(133), quiver(rot90(fliplr(U)),rot90(fliplr(V))), axis('square');axis([1 N 1 N]); drawnow % streamfunction
        end
    end
end