function x = jacobi(FDM, b, tol, x0)

A = FDM.A;
iD = FDM.iD;
R = FDM.R;

x = x0;
err = inf;
while err > tol
    x = iD * (b - R*x);
    res = b - A*x;
    err = norm(res);
end

end