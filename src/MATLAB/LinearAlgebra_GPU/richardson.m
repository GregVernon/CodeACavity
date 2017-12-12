function x = richardson(FDM, b, tol, x0)
weight = 0.2; %2/3;
x = x0;
% f = b - FDM*x;
A = FDM.A;
D = FDM.D; %speye(size(FDM)) .* spdiags(FDM,0);
iD = FDM.iD;

err = inf;
while err > tol
    x = x - weight*(iD*(A*x - b));
    res = b - A*x;
    err = norm(res);
end