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
    [jj,ii] = ind2sub([NY-2,NX-2],node);
    if ii == 1 && jj == 1
        % Bottom Left Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj,ii+1) = R;   %1;
        %         Stencil(jj+1,ii) = T;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = ones(1,3) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj jj+1], [ii ii+1 ii]);
        sVal(sIdx) = [C R T];
    elseif ii == NX-2 && jj == 1
        % Bottom Right Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj,ii-1) = L;   %1;
        %         Stencil(jj+1,ii) = T;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = ones(1,3) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj jj+1], [ii ii-1 ii]);
        sVal(sIdx) = [C L T];
    elseif ii == 1 && jj == NY-2
        % Top Left Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj-1,ii) = B;   %1;
        %         Stencil(jj,ii+1) = R;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = ones(1,3) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj-1 jj], [ii ii ii+1]);
        sVal(sIdx) = [C B R];
    elseif ii == NX-2 && jj == NY-2
        % Top Right Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj,ii-1) = L;   %1;
        %         Stencil(jj-1,ii) = B;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 3];
        sRow(sIdx) = ones(1,3) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj jj-1], [ii ii-1 ii]);
        sVal(sIdx) = [C L B];
    elseif ii == 1
        %             Left Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj,ii+1) = R;   %1;
        %         Stencil(jj+1,ii) = T;   %1;
        %         Stencil(jj-1,ii) = B;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = ones(1,4) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj jj+1 jj-1], [ii ii+1 ii ii]);
        sVal(sIdx) = [C R T B];
    elseif ii == NX-2
        % Right Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj,ii-1) = L;   %1;
        %         Stencil(jj+1,ii) = T;   %1;
        %         Stencil(jj-1,ii) = B;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = ones(1,4) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj jj+1 jj-1], [ii ii-1 ii ii]);
        sVal(sIdx) = [C L T B];
    elseif jj == 1
        % Bottom Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj+1,ii) = T;   %1;
        %         Stencil(jj,ii-1) = L;   %1;
        %         Stencil(jj,ii+1) = R;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = ones(1,4) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj+1 jj jj], [ii ii ii-1 ii+1]);
        sVal(sIdx) = [C T L R];
    elseif jj == NY-2
        % Top Node
        % Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj-1,ii) = B;   %1;
        %         Stencil(jj,ii+1) = R;   %1;
        %         Stencil(jj,ii-1) = L;   %1;
        
        % Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 4];
        sRow(sIdx) = ones(1,4) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj-1 jj jj], [ii ii ii+1 ii-1]);
        sVal(sIdx) = [C B R L];
    else
        % Interior Node
        %         %Slower
        %         Stencil(jj,ii) = C;     %-4;
        %         Stencil(jj,ii+1) = R;   %1;
        %         Stencil(jj,ii-1) = L;   %1;
        %         Stencil(jj+1,ii) = T;   %1;
        %         Stencil(jj-1,ii) = B;   %1;
        
        %Faster
        sIdx = [sIdx(end) + 1 : sIdx(end) + 5];
        sRow(sIdx) = ones(1,5) * node;
        sCol(sIdx) = sub2ind([NY-2 NX-2], [jj jj jj jj+1 jj-1], [ii ii+1 ii-1 ii ii]);
        sVal(sIdx) = [C R L T B];

    end
end

sRow(isnan(sRow)) = [];
sCol(isnan(sCol)) = [];
sVal(isnan(sVal)) = [];
A = sparse(sRow,sCol,sVal,(NY-2)*(NX-2),(NY-2)*(NX-2));
% spy(A)