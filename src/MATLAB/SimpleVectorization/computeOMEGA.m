function OMEGA = computeOMEGA(PSI, OMEGA, dt, h, NY, NX, Visc)

OMEGA0 = OMEGA;
jj=2:NY-1;
ii=2:NX-1;
OMEGA(jj,ii) = OMEGA(jj,ii) + dt * (-0.25*((PSI(jj,ii+1)-PSI(jj,ii-1)).*(OMEGA0(jj+1,ii)-OMEGA0(jj-1,ii)) -(PSI(jj+1,ii)-PSI(jj-1,ii)).*(OMEGA0(jj,ii+1)-OMEGA0(jj,ii-1)))/(h^2) +Visc*(OMEGA0(jj+1,ii)+OMEGA0(jj-1,ii)+OMEGA0(jj,ii+1)+OMEGA0(jj,ii-1)-4.0*OMEGA0(jj,ii))/(h^2)); % vorticity