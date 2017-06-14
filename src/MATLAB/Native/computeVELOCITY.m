function [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, h, NY, NX)
U = zeros(NY,NX);
V = zeros(NY,NX);
VELOCITY = zeros(NY,NX);
cREYNOLDS = zeros(NY,NX);

for jj = 2:NY-1
    for ii = 2:NX-1
        U(jj,ii) = (PSI(jj,ii+1) - PSI(jj,ii-1)) / (2*h);
        V(jj,ii) = (PSI(jj-1,ii) - PSI(jj+1,ii)) / (2*h);
        VELOCITY(jj,ii) = natmat_sqrt(U(jj,ii) * U(jj,ii) + V(jj,ii) * V(jj,ii));
        cREYNOLDS(jj,ii) = natmat_abs(U(jj,ii)) + natmat_abs(V(jj,ii));
    end
end
