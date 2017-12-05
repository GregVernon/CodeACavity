function x = AlternatingAndersonRichardson(FDM, b, tol, x0)
warning ('off','all')
weight = 2/3;
beta = 0.2;
x = x0;

A = FDM.A;
D = FDM.D; %speye(size(FDM)) .* spdiags(FDM,0);
iD = FDM.iD;
I = speye(size(A));


aFreq = 25;
aCount = 0;
nHist = 20;
mCount = 0;
X = zeros(size(x,1),nHist+1);
F = zeros(size(x,1),nHist+1);

iter = 0;
err = inf;
while err > tol
    iter = iter+1;
    aCount = aCount + 1;
    if aCount < aFreq
        if aFreq - (aCount+1) <= nHist
            mCount = mCount+1;
            X(:,mCount) = x;
            F(:,mCount) = f;
        end
        % Richardson extrapolation
        B = weight*I;
        f = iD*(b - A*x);
        x = x + B*f;
        
        if aFreq - (aCount+1) <= nHist
            X(:,mCount) = x - X(:,mCount);
            F(:,mCount) = f - F(:,mCount);
        end
    else
        aCount = 0;
        mCount = 0;
%         B = beta*I - (X + beta*F)*inv(F'*F)*F';
        B = beta*I - (X + beta*F)/(F'*F)*F';
        f = iD*(b - A*x);
        x = x + B*f;    
    end
    
    res = b - A*x;
    err = norm(res);
end
disp(iter)
warning ('on','all')