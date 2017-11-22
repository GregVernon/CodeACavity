function x = weighted_jacobi(FDM, b, tol, x0)

A = FDM.A;
iD = FDM.iD;
R = FDM.R;

weight = 1.1;
x = x0;
err = inf;
iter = 0;
while err > tol
    iter = iter + 1;
    g = iD * (b - R*x);
    f = g - x;
    x = x + weight*f;
    res = b - A*x;
    err = norm(res);
end
% disp(iter)
end