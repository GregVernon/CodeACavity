function PSI = computePSI(PSI, OMEGA, FDM, method, tol, h, NY, NX)
nx = NX;
ny = NY;
jj = 2:NY-1;
ii = 2:NX-1;

b = computeRHS(PSI,OMEGA,NX,NY,h);
x0 = reshape(PSI(jj,ii),(ny-2)*(nx-2),1);
if strcmpi(method,'Direct')
    x = FDM \ b;
elseif strcmpi(method,'CG')
    x = conjgrad(FDM, b, x0, tol);
elseif strcmpi(method,'Decomposition')
    x = FDM \ b;
end

PSI(jj,ii) = reshape(x',ny-2,nx-2);