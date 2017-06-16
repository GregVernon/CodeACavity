function OMEGA = computeOMEGA(PSI, OMEGA, dt, h, NY, NX, Visc)

jj=2:NY-1;
ii=2:NX-1;
OMEGA(jj,ii) = OMEGA(jj,ii) + dt * (-0.25*((PSI(jj,ii+1)-PSI(jj,ii-1)).*(OMEGA(jj+1,ii)-OMEGA(jj-1,ii)) -(PSI(jj+1,ii)-PSI(jj-1,ii)).*(OMEGA(jj,ii+1)-OMEGA(jj,ii-1)))/(h^2) +Visc*(OMEGA(jj+1,ii)+OMEGA(jj-1,ii)+OMEGA(jj,ii+1)+OMEGA(jj,ii-1)-4.0*OMEGA(jj,ii))/(h^2)); % vorticity