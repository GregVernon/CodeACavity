function x = jacobi_Mapped(FDM, b, tol, x0, dx, dy, NX, NY)

A = FDM.A;
% iD = FDM.iD;
iD = 1/(-2*(dy^2 + dx^2));
R = FDM.R;

x = x0;
err = inf;

bV = b;
bM = reshape(b,NY,NX);
x = reshape(x,NY,NX);
jj = 2:NY-1;
ii = 2:NX-1;
while err > tol
%     x = iD * (b - R*x);
%     Rx = ((dx^2*x(jj+1,ii)) + (dx^2*x(jj-1,ii)) + (dy^2*x(jj,ii+1)) + (dy^2*x(jj,ii-1)));
%     x(jj,ii) =  iD * (bM(jj,ii) - Rx);
    x(jj,ii) = arrayfun(@xArrayFun, x(jj+1,ii), x(jj-1,ii), x(jj,ii+1), x(jj,ii-1), bM(jj,ii), iD, dx, dy);
    res = bV - A*x(:);
    err = norm(res);
end

x = reshape(x,NY*NX,1);
end

function xji = xArrayFun(xjp1i, xjm1i, xjip1, xjim1, bji, iD, dx, dy)
    
    Rx = ((dx^2*xjp1i) + (dx^2*xjm1i) + (dy^2*xjip1) + (dy^2*xjim1));
    xji = iD * (bji - Rx);
end