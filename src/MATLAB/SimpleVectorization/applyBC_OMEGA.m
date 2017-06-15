function OMEGA = applyBC_OMEGA(PSI, OMEGA, h, NY, NX)

for jj = 2:NY-1
    OMEGA(jj,1) = -2.0*PSI(jj,2)/(h^2); % vorticity on bottom wall
    OMEGA(jj,NX)= -2.0*PSI(jj,NX-1)/(h^2) - 2.0/h; % vorticity on top wall
end

for ii = 2:NX-1
    OMEGA(1,ii) = -2.0*PSI(2,ii)/(h^2);       % vorticity on right wall
    OMEGA(NY,ii)= -2.0*PSI(NY-1,ii)/(h^2); % vorticity on left wall
end