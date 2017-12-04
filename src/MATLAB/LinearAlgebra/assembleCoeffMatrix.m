function A = assembleCoeffMatrix(dx, dy, NX, NY)
%% Generate A for a domain from 1:N
sRow = nan(1, 5 * (NY-2) * (NX-2));
sCol = nan(1, 5 * (NY-2) * (NX-2));
sVal = nan(1, 5 * (NY-2) * (NX-2));
sIdx = 0;

% syms C L R D U
C = -2*(dy^2 + dx^2); %-4;
L = 1 * dy^2; %1;
R = 1 * dy^2; %1;
D = 1 * dx^2; %1;
U = 1 * dx^2; %1;
nNodes = (NX-2)*(NY-2);
for node = 1:nNodes
    jj = rem(node-1, (NY-2)) + 1;
    ii = (node - jj)/(NY-2) + 1;
    if ii == 1 && jj == 1
        % Bottom Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj jj+1] + ([ii ii+1 ii] - 1).*(NY-2);
        sVal(sIdx) = [C R U];
    elseif ii == NX-2 && jj == 1
        % Bottom Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj jj+1] + ([ii ii-1 ii] - 1).*(NY-2);
        sVal(sIdx) = [C L U];
    elseif ii == 1 && jj == NY-2
        % Top Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj-1 jj] + ([ii ii ii+1] - 1).*(NY-2);
        sVal(sIdx) = [C D R];
    elseif ii == NX-2 && jj == NY-2
        % Top Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj jj-1] + ([ii ii-1 ii] - 1).*(NY-2);
        sVal(sIdx) = [C L D];
    elseif ii == 1
        % Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj jj+1 jj-1] + ([ii ii+1 ii ii] - 1).*(NY-2);
        sVal(sIdx) = [C R U D];
    elseif ii == NX-2
        % Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj jj+1 jj-1] + ([ii ii-1 ii ii] - 1).*(NY-2);
        sVal(sIdx) = [C L U D];
    elseif jj == 1
        % Bottom Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj+1 jj jj] + ([ii ii ii-1 ii+1] - 1).*(NY-2);
        sVal(sIdx) = [C U L R];
    elseif jj == NY-2
        % Top Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj-1 jj jj] + ([ii ii ii+1 ii-1] - 1).*(NY-2);
        sVal(sIdx) = [C D R L];
    else
        % Interior Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 5];
        sRow(sIdx) = [node node node node node];
        sCol(sIdx) = [jj jj jj jj+1 jj-1] + ([ii ii+1 ii-1 ii ii] - 1).*(NY-2); 
        sVal(sIdx) = [C R L U D];
    end
end

sRow(isnan(sRow)) = [];
sCol(isnan(sCol)) = [];
sVal(isnan(sVal)) = [];
A = sparse(sRow,sCol,sVal,(NY-2)*(NX-2),(NY-2)*(NX-2));
% spy(A)