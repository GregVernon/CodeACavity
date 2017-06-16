function [PSI,adiffPSI] = jacobiIteration(PSI_jp1,PSI_jm1,PSI_ip1,PSI_im1,PSI0,OMEGA,h)

PSI = (PSI_jp1 + PSI_jm1 + PSI_ip1 + PSI_im1)./4 + (h^2.*OMEGA);
adiffPSI = PSI - PSI0;