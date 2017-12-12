function [x] = conjgrad(A, b, x, tol)
    r = b - A * x;
    p = r;
    rsold = r' * r;
    res = inf;
    while res >= tol
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        res = r' * r;
        
        p = r + (res / rsold) * p;
        rsold = res;
    end
end