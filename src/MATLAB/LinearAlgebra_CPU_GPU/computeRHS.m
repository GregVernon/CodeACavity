function b = computeRHS(PSI_jp1,PSI_jm1,PSI_ip1,PSI_im1,OMEGA,NX,NY,h,ii,jj)
%% Generate b for a domain from 2:N-1
b = 0;
%% Bottom Left Node
if jj == 2 && ii == 2
    b = -(h^2 *  OMEGA) + (PSI_im1 + PSI_jm1);
elseif jj == 2 && ii == NX-1
    % Bottom Right Node
    b = -(h^2 *  OMEGA) + (PSI_ip1 - PSI_jm1);
elseif jj == NY-1 && ii == 2
    % Top Left Node
    b = -(h^2 *  OMEGA) - (PSI_jp1 + PSI_im1);
elseif jj == NY-1 && ii == NX-1
    % Top Right Node
    b = -(h^2 *  OMEGA) - (PSI_jp1 + PSI_ip1);
elseif jj > 2 && jj < NY-1 && ii == 2
    % Left Node
    b = -(h^2 *  OMEGA) - (PSI_im1);
elseif jj > 2 && jj < NY-1 && ii == NX-1
    % Right Node
    b = -(h^2 *  OMEGA) - (PSI_ip1);
elseif jj == 2 && ii > 2 && ii < NX-1
    % Bottom Node
    b = -(h^2 *  OMEGA) - (PSI_jm1);
elseif jj == NY-1 && ii > 2 && ii < NX-1
    % Top Node
    b = -(h^2 *  OMEGA) - (PSI_jp1);
elseif jj >= 3 && jj <= NY-2 && ii >= 3 && ii <= NX-2
    % Interior Node
    b = -(h^2 * OMEGA);
end
% Reshape into a vector
% b = reshape(b(2:NY-1,2:NX-1),(NX-2)*(NY-2),1);