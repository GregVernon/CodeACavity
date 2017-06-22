function [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI_jp1,PSI_jm1,PSI_ip1,PSI_im1, h)
U = (PSI_ip1 - PSI_im1) / (2*h);
V = (PSI_jm1 - PSI_jp1) / (2*h);
VELOCITY = sqrt(U^2 + V^2);
cREYNOLDS = abs(U) + abs(V);