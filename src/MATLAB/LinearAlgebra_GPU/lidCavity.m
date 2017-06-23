function lidCavity(N, tmax, tol, plotFlow)

NY = N;
NX = N;
jj = ((2:NY-1)');
ii = ((2:NX-1));
h=1.0/(NY-1);

Visc=.01;

dt_max = tmax / 1000;

PSI = zeros(NY,NX,'gpuArray');
OMEGA = zeros(NY,NX,'gpuArray');
OMEGA(jj,[NX-1 NX])= -2.0/h;
U = zeros(NY,NX,'gpuArray');
V = zeros(NY,NX,'gpuArray');

VELOCITY = zeros(NY,NX,'gpuArray');
cREYNOLDS = zeros(NY,NX,'gpuArray');


FDM = assembleCoeffMatrix(NX,NY);
method = 'cg';
t=0.0;
pIter = 0;
while t < tmax % start the time integration
    pIter = pIter+1;
    % Compute timestep
    dt = delta_t(VELOCITY, cREYNOLDS, dt_max, h);
    % Compute streamfunction
    PSI = computePSI(PSI, OMEGA, FDM, method,tol, h, NY, NX);
    % Apply vorticity boundary conditions
    OMEGA = applyBC_OMEGA(PSI, OMEGA, h, NY, NX);
    % Compute vorticity
    OMEGA(jj,ii) = arrayfun(@computeOMEGA,PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1), OMEGA(jj+1,ii), OMEGA(jj-1,ii), OMEGA(jj,ii+1), OMEGA(jj,ii-1), OMEGA(jj,ii), dt, h, Visc);
    % Compute Velocity
    [U(jj,ii), V(jj,ii), VELOCITY(jj,ii), cREYNOLDS(jj,ii)] = arrayfun(@computeVELOCITY,PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1), h);
    % Increment time value by timestep
    t=t+dt;
    %% plot
    if pIter == 1e2
        pIter = 0;
        disp(['Time: ' num2str(t)])
        if plotFlow == true
            subplot(131), contourf(rot90(fliplr(OMEGA))), axis('square'); colorbar% plot vorticity
            subplot(132), contourf(rot90(fliplr(PSI))), axis('square'); colorbar% streamfunction
            subplot(133), quiver(rot90(fliplr(U)),rot90(fliplr(V))), axis('square');axis([1 N 1 N]); drawnow % streamfunction
        end
    end
end