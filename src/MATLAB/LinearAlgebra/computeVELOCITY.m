function [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, dx, dy, NY, NX)
U = zeros(NY,NX);
V = zeros(NY,NX);



jj = 2:NY-1;
ii = 2:NX-1;
U(jj,ii) =  (PSI(jj+1,ii) - PSI(jj-1,ii)) / (2*dy);
V(jj,ii) = -(PSI(jj,ii+1) - PSI(jj,ii-1)) / (2*dx);
VELOCITY = sqrt(U.^2 + V.^2);
cREYNOLDS = abs(U) + abs(V);