function [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, h, NY, NX)
U = zeros(NY,NX);
V = zeros(NY,NX);
VELOCITY = zeros(NY,NX);
cREYNOLDS = zeros(NY,NX);

for jj = 2:NY-1
    for ii = 2:NX-1
        U(jj,ii) = (PSI(jj,ii+1) - PSI(jj,ii-1)) / (2*h);
        V(jj,ii) = (PSI(jj-1,ii) - PSI(jj+1,ii)) / (2*h);
        VELOCITY(jj,ii) = sqrt(U(jj,ii)^2 + V(jj,ii)^2);
        cREYNOLDS(jj,ii) = abs(U(jj,ii)) + abs(V(jj,ii));
    end
end
