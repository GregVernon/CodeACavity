function PSI = computePSI(PSI, OMEGA, h, NY, NX, tol)
jj = gpuArray(2:NY-1);
ii = gpuArray(2:NX-1);
err = gpuArray(inf);
kPSI = 0;
while err > tol
    kPSI = kPSI + 1;
    PSI0 = PSI;
    PSI(jj,ii)=0.25*(PSI(jj+1,ii)+PSI(jj-1,ii)+PSI(jj,ii+1)+PSI(jj,ii-1)+h^2*OMEGA(jj,ii));
    err = sum(sum(abs(PSI0 - PSI)));
end