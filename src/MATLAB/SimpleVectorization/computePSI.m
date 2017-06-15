function PSI = computePSI(PSI, OMEGA, h, NY, NX, tol)
err = inf;
while err > tol
    err = 0;
    PSI0=PSI;
    for jj=2:NY-1
        for ii=2:NX-1
            PSI(jj,ii)=0.25*(PSI(jj+1,ii)+PSI(jj-1,ii)+PSI(jj,ii+1)+PSI(jj,ii-1)+h^2*OMEGA(jj,ii));
            err = err + abs(PSI0(jj,ii) - PSI(jj,ii));
        end
    end
end