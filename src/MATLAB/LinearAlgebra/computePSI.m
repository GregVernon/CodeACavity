function PSI = computePSI(PSI, OMEGA, FDM, method, tol, dx, dy, NX, NY)

ii = 1:NX;
jj = 1:NY;

b = computeRHS(PSI,OMEGA,NX,NY,dx,dy);
x0 = reshape(PSI(jj,ii),(NX)*(NY),1);

if strcmpi(method,'Direct')
    x = FDM \ b;
elseif strcmpi(method,'CG')
    x = conjgrad(FDM, b, x0, tol);
elseif strcmpi(method,'Decomposition')
    x = FDM \ b;
elseif strcmpi(method,'Jacobi')
    x = jacobi(FDM,b,tol,x0);
elseif strcmpi(method,'Mapped Jacobi')
    x = jacobi_Mapped(FDM,b,tol,x0,dx,dy,NX,NY);
elseif strcmpi(method,'Weighted Jacobi')
    x = weighted_jacobi(FDM,b,tol,x0);
elseif strcmpi(method,'Richardson')
    x = richardson(FDM,b,tol,x0);
elseif strcmpi(method,'Mapped Richardson')
    x = richardson_Mapped(FDM,b,tol,x0,dx,dy,NX,NY);
elseif strcmpi(method,'AlternatingAndersonRichardson')
    x = AlternatingAndersonRichardson(FDM,b,tol,x0);
end

PSI(jj,ii) = reshape(x',NY,NX);