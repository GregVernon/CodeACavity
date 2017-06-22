function PSI = computePSI(PSI, OMEGA, FDM, method,tol, h, NY, NX)
nx = NX;
ny = NY;
jj = [2:NY-1]';
ii = 2:NX-1;

b = arrayfun(@computeRHS,PSI(jj+1,ii),PSI(jj-1,ii),PSI(jj,ii+1),PSI(jj,ii-1),OMEGA(jj,ii),NX,NY,h,ii,jj);
b = reshape(b,(NY-2)*(NX-2),1);
b = gather(b);
x0 = reshape(PSI(jj,ii),(ny-2)*(nx-2),1);
x0 = gather(x0);

if strcmpi(method,'Direct')
    x = FDM\b;
elseif strcmpi(method,'CG')
    x = conjgrad(FDM, b, x0, tol);
end

x = gpuArray(x);
PSI(jj,ii) = reshape(x',ny-2,nx-2);