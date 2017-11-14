function A = assembleCoeffMatrix(NX, NY)
%% Generate A for a domain from 1:N
sRow = nan(1, 5 * (NY-2) * (NX-2));
sCol = nan(1, 5 * (NY-2) * (NX-2));
sVal = nan(1, 5 * (NY-2) * (NX-2));
sIdx = 0;

% syms C L R B T
C = -4;%*(dy^2 + dx^2); %-4;
L = 1; %1;
R = 1; %1;
B = 1; %1;
T = 1; %1;
nNodes = (NX-2)*(NY-2);
for node = 1:nNodes
    jj = rem(node-1, (NY-2)) + 1;
    ii = (node - jj)/(NY-2) + 1;
    if ii == 1 && jj == 1
        % Bottom Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj jj+1] + ([ii ii+1 ii] - 1).*(NY-2);
        sVal(sIdx) = [C R T];
    elseif ii == NX-2 && jj == 1
        % Bottom Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj jj+1] + ([ii ii-1 ii] - 1).*(NY-2);
        sVal(sIdx) = [C L T];
    elseif ii == 1 && jj == NY-2
        % Top Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj-1 jj] + ([ii ii ii+1] - 1).*(NY-2);
        sVal(sIdx) = [C B R];
    elseif ii == NX-2 && jj == NY-2
        % Top Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = [node node node];
        sCol(sIdx) = [jj jj jj-1] + ([ii ii-1 ii] - 1).*(NY-2);
        sVal(sIdx) = [C L B];
    elseif ii == 1
        % Left Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj jj+1 jj-1] + ([ii ii+1 ii ii] - 1).*(NY-2);
        sVal(sIdx) = [C R T B];
    elseif ii == NX-2
        % Right Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj jj+1 jj-1] + ([ii ii-1 ii ii] - 1).*(NY-2);
        sVal(sIdx) = [C L T B];
    elseif jj == 1
        % Bottom Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj+1 jj jj] + ([ii ii ii-1 ii+1] - 1).*(NY-2);
        sVal(sIdx) = [C T L R];
    elseif jj == NY-2
        % Top Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = [node node node node];
        sCol(sIdx) = [jj jj-1 jj jj] + ([ii ii ii+1 ii-1] - 1).*(NY-2);
        sVal(sIdx) = [C B R L];
    else
        % Interior Node
        sIdx = [sIdx(end) + 1 : sIdx(end) + 5];
        sRow(sIdx) = [node node node node node];
        sCol(sIdx) = [jj jj jj jj+1 jj-1] + ([ii ii+1 ii-1 ii ii] - 1).*(NY-2); 
        sVal(sIdx) = [C R L T B];
    end
end

sRow(isnan(sRow)) = [];
sCol(isnan(sCol)) = [];
sVal(isnan(sVal)) = [];
A = sparse(sRow,sCol,sVal,(NY-2)*(NX-2),(NY-2)*(NX-2));
% spy(A)