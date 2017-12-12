function x = richardson_Mapped(FDM, b, tol, x0,dx,dy,NX,NY)
weight = 2/3;
x = x0;
% f = b - FDM*x;
A = FDM.A;
D = FDM.D; %speye(size(FDM)) .* spdiags(FDM,0);
% iD = FDM.iD;
D = (-2*(dy^2 + dx^2));
iD = 1/D;

err = inf;
bV = b;
bM = reshape(b,NY,NX);
x = reshape(x,NY,NX);
jj = 2:NY-1;
ii = 2:NX-1;
while err > tol
%     Ax = ((dx^2*x(jj+1,ii)) + (dx^2*x(jj-1,ii)) + (dy^2*x(jj,ii+1)) + (dy^2*x(jj,ii-1))) + (-2*(dy^2 + dx^2)*x(jj,ii));
%     x(jj,ii) = x(jj,ii) + weight*(iD*(b(jj,ii) - Ax));
    x(jj,ii) = arrayfun(@xArrayFun, x(jj,ii), x(jj+1,ii), x(jj-1,ii), x(jj,ii+1), x(jj,ii-1), bM(jj,ii), iD, weight, dx, dy);
    res = bV - A*x(:);
    err = norm(res);
end

x = reshape(x,NY*NX,1);
end

function xji = xArrayFun(xji, xjp1i, xjm1i, xjip1, xjim1, bji, iD, weight, dx, dy)
    Ax = ((dx^2*xjp1i) + (dx^2*xjm1i) + (dy^2*xjip1) + (dy^2*xjim1)) + (-2*(dy^2 + dx^2)*xji);
    xji = xji + weight*(iD*(bji - Ax));
end