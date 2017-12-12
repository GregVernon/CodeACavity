function OMEGA = applyBC_OMEGA(PSI, OMEGA, dx, dy, NY, NX)

ii = 2:NX-1;
OMEGA(1,ii) = -2.0*PSI(2,ii)/(dy^2); % vorticity on bottom wall
OMEGA(NY,ii)= -2.0*PSI(NY-1,ii)/(dy^2) - 2.0/dy; % vorticity on top wall

jj = 2:NY-1;
OMEGA(jj,1) = -2.0*PSI(jj,2)/(dx^2); % vorticity on left wall
OMEGA(jj,NX)= -2.0*PSI(jj,NX-1)/(dx^2); % vorticity on right wall