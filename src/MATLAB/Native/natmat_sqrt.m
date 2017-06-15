function x = natmat_sqrt(y0)

tol = 1e-14;
if y0 == 0
    x = 0;
    return
elseif y0 < 1
    x = 1;
else
    x = y0;
end

err = x * x;
while err >= tol
    x = x - (x * x - y0) / (2 * x);
    err = natmat_abs((x * x) - y0) / y0;
end