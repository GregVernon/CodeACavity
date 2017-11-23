function A = assembleCoeffMatrix(NX, NY)
%% Generate A for a domain from 1:N
sRow = nan(1, 5 * NY * NX);
sCol = nan(1, 5 * NY * NX);
sVal = nan(1, 5 * NY * NX);
sIdx = 0;

% syms C L R B T
C = -4;%*(dy^2 + dx^2); %-4;
L = 1; %1;
R = 1; %1;
B = 1; %1;
T = 1; %1;
nNodes = NX * NY;
for node = 1:nNodes
    jj = rem(node-1, NY) + 1;
    ii = (node - jj)/(NY) + 1;
    if ii == 1 && jj == 1
        % Bottom Left Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + ((ii - 1) .* NY);
        sVal(sIdx) = 1;
    elseif ii == NX && jj == 1
        % Bottom Right Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        ssCol(sIdx) = jj + ((ii - 1) .* NY);
        sVal(sIdx) = 1;
    elseif ii == 1 && jj == NY
        % Top Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        ssCol(sIdx) = [jj jj+1 jj+2] + (([ii ii ii] - 1) .* NY);
        sVal(sIdx) = [-3 4 -1];
    elseif ii == NX && jj == NY
        % Top Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        ssCol(sIdx) = [jj jj+1 jj+2] + (([ii ii ii] - 1) .* NY);
        sVal(sIdx) = [-3 4 -1];
    elseif ii == 1
        % Left Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + ((ii - 1) .* NY);
        sVal(sIdx) = 1;
    elseif ii == NX-2
        % Right Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + ((ii - 1) .* NY);
        sVal(sIdx) = 1;
    elseif jj == 1
        % Bottom Node
        sIdx = sIdx(end) + 1;
        sRow(sIdx) = node;
        sCol(sIdx) = jj + ((ii - 1) .* NY);
        sVal(sIdx) = 1;
    elseif jj == NY-2
        % Top Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        ssCol(sIdx) = [jj jj+1 jj+2] + (([ii ii ii] - 1) .* NY);
        sVal(sIdx) = [-3 4 -1];
    else
        % Interior Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 5];
        sRow(sIdx) = [node node node node node];
        sCol(sIdx) = [jj jj jj jj+1 jj-1] + (([ii ii+1 ii-1 ii ii] - 1).* NY);
        sVal(sIdx) = [C R L T B];
    end
end

sRow(isnan(sRow)) = [];
sCol(isnan(sCol)) = [];
sVal(isnan(sVal)) = [];
A = sparse(sRow,sCol,sVal,NX*NY,NX*NY);
% spy(A)
