function [U, V, VELOCITY, cREYNOLDS] = computeVELOCITY(PSI, dx, dy, NY, NX)
U = zeros(NY,NX,'gpuArray');
V = zeros(NY,NX,'gpuArray');



jj = 2:NY-1;
ii = 2:NX-1;
ji = jj + (ii - 1).*(NY);

jj = [2:NY-1]+1;
ii = [2:NX-1];
jp1i = jj + (ii - 1).*(NY);

jj = [2:NY-1]-1;
ii = [2:NX-1];
jm1i = jj + (ii - 1).*(NY);

jj = [2:NY-1];
ii = [2:NX-1]+1;
jip1 = jj + (ii - 1).*(NY);

jj = [2:NY-1];
ii = [2:NX-1]-1;
jim1 = jj + (ii - 1).*(NY);

% U(jj,ii) =  arrayfun(@uArrayFun, PSI(jj+1,ii), PSI(jj-1,ii), dy); %(PSI(jj+1,ii) - PSI(jj-1,ii)) / (2*dy);
% V(jj,ii) =  arrayfun(@vArrayFun, PSI(jj,ii+1), PSI(jj,ii-1), dx); %-(PSI(jj,ii+1) - PSI(jj,ii-1)) / (2*dx);
U(ji) =   (PSI(jp1i) - PSI(jm1i)) / (2*dy); % arrayfun(@uArrayFun, PSI(jp1i), PSI(jm1i), dy);
V(ji) =  -(PSI(jip1) - PSI(jim1)) / (2*dx); % arrayfun(@vArrayFun, PSI(jip1), PSI(jim1), dx);
VELOCITY =  sqrt(U.^2 + V.^2);
cREYNOLDS = abs(U) + abs(V);

end

function Uji = uArrayFun(PSIjp1i, PSIjm1i, dy)
    Uji = (PSIjp1i - PSIjm1i) / (2*dy);
end

function Vji = vArrayFun(PSIjip1, PSIjim1, dx)
    Vji = -(PSIjip1 - PSIjim1) / (2*dx);
end

