function A = assembleCoeffMatrix(dx, dy, NX, NY)
%% Generate A for a domain from 1:N
sRow = nan(1, 5 * (NY) * (NX));
sCol = nan(1, 5 * (NY) * (NX));
sVal = nan(1, 5 * (NY) * (NX));
sIdx = 0;

% syms C L R D U
I = 1;
C = -2*(dy^2 + dx^2); %-4;
L = 1 * dy^2; %1;
R = 1 * dy^2; %1;
D = 1 * dx^2; %1;
U = 1 * dx^2; %1;
nNodes = (NX)*(NY);
for node = 1:nNodes
    jj = rem(node-1, NY) + 1;
    ii = (node - jj)/(NY) + 1;
    if ii == 1 && jj == 1
        % Bottom Left Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif ii == NX && jj == 1
        % Bottom Right Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif ii == 1 && jj == NY
        % Top Left Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif ii == NX && jj == NY
        % Top Right Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif ii == 1
        % Left Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif ii == NX
        % Right Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif jj == 1
        % Bottom Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    elseif jj == NY
        % Top Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + (ii - 1).*(NY);
        sVal(sIdx) = I;
    else
        % Interior Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 5];
        sRow(sIdx) = [node node node node node];
        sCol(sIdx) = [jj jj jj jj+1 jj-1] + ([ii ii+1 ii-1 ii ii] - 1).*(NY); 
        sVal(sIdx) = [C R L U D];
    end
end

sRow(isnan(sRow)) = [];
sCol(isnan(sCol)) = [];
sVal(isnan(sVal)) = [];
A = sparse(sRow,sCol,sVal,(NY)*(NX),(NY)*(NX));
% spy(A)