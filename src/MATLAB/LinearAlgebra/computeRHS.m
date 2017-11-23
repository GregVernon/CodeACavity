function b = computeRHS(PSI,OMEGA,U,NX,NY,h)
dx = h;
dy = h;
%% Preallocate b
b = zeros(NX*NY,1);
%% Bottom Left Node
jj = 1;
ii = 1;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = 0;
%% Bottom Right Node
jj = 1;
ii = NX;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = 0;
%% Top Left Node
jj = NY;
ii = 1;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = 2*dy * U(jj,ii);
%% Top Right Node
jj = NY;
ii = NX;
bidx = sub2ind([NY,NX],jj,ii);
b(bidx) = 2*dy * U(jj,ii);
%% Left Node
jj = 1:NY;
ii = 1;
bidx = sub2ind([NY,NX],jj,ii*ones(size(jj)));
b(bidx) = 0;
%% Right Node
jj = 1:NY;
ii = NX;
bidx = sub2ind([NY,NX],jj,ii*ones(size(jj)));
b(bidx) = 0;
%% Bottom Node
jj = 1;
ii = 1:NX;
bidx = sub2ind([NY,NX],jj*ones(size(ii)),ii);
b(bidx) = 0;
%% Top Node
jj = NY;
ii = 1:NX;
bidx = sub2ind([NY,NX],jj*ones(size(ii)),ii);
b(bidx) = 2*dy * U(jj,ii);
%% Interior Node
jj = 2:NY-1;
ii = 2:NX-1;
JJ = reshape(repmat(jj',1,NX-2),(NY-2)*(NX-2),1);
II = reshape(repmat(ii, NY-2,1),(NY-2)*(NX-2),1);
bidx = sub2ind([NY,NX],JJ,II);
b(bidx) = -(dx^2*dy^2 * OMEGA(jj,ii));
