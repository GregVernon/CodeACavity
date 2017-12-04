function PSI = computePSI(PSI, OMEGA, FDM, method, tol, dx, dy, NX, NY)
ii = 2:NX-1;
jj = 2:NY-1;

b = computeRHS(PSI,OMEGA,NX,NY,dx,dy);
x0 = reshape(PSI(jj,ii),(NX-2)*(NY-2),1);
if strcmpi(method,'Direct')
    x = FDM \ b;
elseif strcmpi(method,'CG')
    x = conjgrad(FDM, b, x0, tol);
elseif strcmpi(method,'Decomposition')
    x = FDM \ b;
elseif strcmpi(method,'Jacobi')
    x = jacobi(FDM,b,tol,x0);
end

PSI(jj,ii) = reshape(x',NY-2,NX-2);