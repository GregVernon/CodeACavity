function b = computeRHS(PSI,OMEGA,NX,NY,dx,dy)
%% Generate b for a domain from 2:N-1
b = zeros(NX*NY,1);
%% Bottom Left Node
jj = 2;
ii = 2;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dx^2 * PSI(jj-1,ii)) - (dy^2 * PSI(jj,ii-1));
%% Bottom Right Node
jj = 2;
ii = NX-1;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dx^2 * PSI(jj-1,ii)) - (dy^2 * PSI(jj,ii+1));
%% Top Left Node
jj = NY-1;
ii = 2;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dx^2 * PSI(jj+1,ii)) - (dy^2 * PSI(jj,ii-1));
%% Top Right Node
jj = NY-1; 
ii = NX-1;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dx^2 * PSI(jj+1,ii)) - (dy^2 * PSI(jj,ii+1));
%% Left Node
jj = 2:NY-1;
ii = 2;
bidx = sub2ind([NY,NX],jj,ii*ones(size(jj)));
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dy^2 * PSI(jj,ii-1));
%% Right Node
jj = 2:NY-1;
ii = NX-1;
bidx = sub2ind([NY,NX],jj,ii*ones(size(jj)));
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dy^2 * PSI(jj,ii+1));
%% Bottom Node
jj = 2;
ii = 2:NX-1; 
bidx = sub2ind([NY,NX],jj*ones(size(ii)),ii);
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dx^2 * PSI(jj-1,ii));
%% Top Node
jj = NY-1;
ii = 2:NX-1;
bidx = sub2ind([NY,NX],jj*ones(size(ii)),ii);
b(bidx) = -(dx^2*dy^2 *  OMEGA(jj,ii)) - (dx^2 * PSI(jj+1,ii));
%% Interior Node
jj = 3:NY-2;
ii = 3:NX-2;
JJ = reshape(repmat([3:NY-2]',1,NX-4),(NY-4)*(NX-4),1);
II = reshape(repmat(3:NX-2,NY-4,1),(NY-4)*(NX-4),1);
bidx = sub2ind([NY,NX],JJ,II);
b(bidx) = -(dx^2*dy^2 * OMEGA(jj,ii));
%% Remove Boundary Nodes from RHS vector
jj = [[1:NY] [1:NY] [1*ones(1,NX)] [NY*ones(1,NX)]];
ii = [[1*ones(1,NY)] [NX*ones(1,NY)] [1:NX] [1:NX]];
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = [];
