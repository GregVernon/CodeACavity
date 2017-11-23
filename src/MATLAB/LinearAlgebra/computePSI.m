function PSI = computePSI(PSI, OMEGA, FDM, method, tol, h, NY, NX)

b = computeRHS(PSI,OMEGA,NX,NY,h);
x0 = reshape(PSI,NX*NY,1);
if strcmpi(method,'Direct')
    x = FDM \ b;
elseif strcmpi(method,'CG')
    x = conjgrad(FDM, b, x0, tol);
elseif strcmpi(method,'Decomposition')
    x = FDM \ b;
elseif strcmpi(method,'Jacobi')
    x = jacobi(FDM,b,tol,x0);
end

PSI = reshape(x',NY,NX);
