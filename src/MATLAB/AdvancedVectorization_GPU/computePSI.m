function PSI = computePSI(PSI, OMEGA, h, NY, NX, tol)
jj = gpuArray(2:NY-1);
ii = gpuArray(2:NX-1);
err = gpuArray(inf);
kPSI = 0;
while err > tol
    kPSI = kPSI + 1;
    [PSI(jj,ii), adiffPSI] = arrayfun(@jacobiIteration,PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),PSI(jj,ii),OMEGA(jj,ii),h);
    err = sum(sum(adiffPSI));
end