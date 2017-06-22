function OMEGA = applyBC_OMEGA(PSI_jp1,PSI_jm1,PSI_ip1,PSI_im1, OMEGA, h,JJ,II, NY, NX)

if II == 1 && (JJ >= 2 && JJ <= NY-1)
    OMEGA = -2.0*PSI_ip1/(h^2); % vorticity on bottom wall
end

if II == NX && (JJ >= 2 && JJ <= NY-1)
    OMEGA= -2.0*PSI_im1/(h^2) - 2.0/h; % vorticity on top wall
end

if JJ == 1 && (II >= 2 && II <= NX-1)
    OMEGA = -2.0*PSI_jp1/(h^2);       % vorticity on right wall
end

if JJ == NY && (II >= 2 && II <= NX-1)
    OMEGA = -2.0*PSI_jm1/(h^2); % vorticity on left wall
end
