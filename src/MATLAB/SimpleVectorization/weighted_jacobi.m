function PSI = weighted_jacobi(PSI, OMEGA, h, NY, NX, tol)
jj = 2:NY-1;
ii = 2:NX-1;
weight = 2/3; % Typical
err = inf;
while err > tol
    PSI0=PSI;
    PSI(jj,ii)= (1 - weight)* PSI(jj,ii) + weight*(0.25*(PSI(jj+1,ii)+PSI(jj-1,ii)+PSI(jj,ii+1)+PSI(jj,ii-1)+h^2*OMEGA(jj,ii)));
    err = sum(sum(abs(PSI0 - PSI)));
end