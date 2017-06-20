function [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, h, NY, NX)
U = gpuArray(zeros(NY,NX));
V = gpuArray(zeros(NY,NX));



jj = 2:NY-1;
ii = 2:NX-1;
U(jj,ii) = (PSI(jj,ii+1) - PSI(jj,ii-1)) / (2*h);
V(jj,ii) = (PSI(jj-1,ii) - PSI(jj+1,ii)) / (2*h);
VELOCITY = bsxfun(@hypot,U,V);%sqrt(U.^2 + V.^2);
cREYNOLDS = abs(U) + abs(V);