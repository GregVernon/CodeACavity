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
b = reshape(b,NY,NX);
x = reshape(x,NY,NX);
jj = 2:NY-1;
ii = 2:NX-1;
while err > tol
    x0 = x;
    Ax = ((dx^2*x(jj+1,ii)) + (dx^2*x(jj-1,ii)) + (dy^2*x(jj,ii+1)) + (dy^2*x(jj,ii-1))) + (-2*(dy^2 + dx^2)*x(jj,ii));
    x(jj,ii) = x(jj,ii) + weight*(iD*(b(jj,ii) - Ax));
%     res = b - A*x;
    err = sum(sum(abs(x-x0)));
end

x = reshape(x,NY*NX,1);
