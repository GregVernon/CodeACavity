function x = weighted_jacobi(FDM, b, tol, x0)

A = FDM.A;
iD = FDM.iD;
R = FDM.R;

weight = 2/3;
x = x0;
err = inf;
iter = 0;
while err > tol
    iter = iter + 1;
    x = weight*iD * (b - R*x) + (1-weight)*x;
    res = b - A*x;
    err = norm(res);
end
% disp(iter)
end