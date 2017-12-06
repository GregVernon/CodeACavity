function x = jacobi_Mapped(FDM, b, tol, x0, dx, dy, NX, NY)

A = FDM.A;
% iD = FDM.iD;
iD = 1/(-2*(dy^2 + dx^2));
R = FDM.R;

x = x0;
err = inf;

b = reshape(b,NY,NX);
x = reshape(x,NY,NX);
jj = 2:NY-1;
ii = 2:NX-1;
while err > tol
%     x = iD * (b - R*x);
    x0 = x;
    Rx = ((dx^2*x(jj+1,ii)) + (dx^2*x(jj-1,ii)) + (dy^2*x(jj,ii+1)) + (dy^2*x(jj,ii-1)));
    x(jj,ii) =  iD * (b(jj,ii) - Rx);
%     res = b - A*x;
%     res = x - x0
    err = sum(sum(abs(x-x0)));
end

x = reshape(x,NY*NX,1);
end