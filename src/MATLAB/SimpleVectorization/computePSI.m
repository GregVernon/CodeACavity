function PSI = computePSI(PSI, OMEGA, h, NY, NX, method, tol)

if strcmpi(method,'Jacobi')
    PSI = jacobi(PSI,OMEGA,h,NY,NX,tol);
elseif strcmpi(method, 'Weighted-Jacobi')
    PSI = weighted_jacobi(PSI,OMEGA,h,NY,NX,tol);
elseif strcmpi(method, 'SR-Jacobi')
elseif strcmpi(method, 'Gauss-Seidel')
elseif strcmpi(method, 'SOR-Gauss-Seidel')
end

