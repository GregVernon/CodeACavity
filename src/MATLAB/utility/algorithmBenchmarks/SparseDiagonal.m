clear
maxNumCompThreads(1);
close all
warning('off','all')
n = 2:8;
N = 2.^n;
for ii = 1:length(N)
    disp(n(ii))
    A = gallery('poisson',N(ii));
    matSize(ii) = size(A,1);
    dVal = spdiags(A,0);
    I = speye(size(A));
    [iRow,iCol] = find(I);
    [B,~] = spdiags(A);
    D1 = @() sparse(iRow,iCol,dVal,size(A,1),size(A,2));
    D2 = @() spdiags(B(:,3),0,size(A,1),size(A,2));
    D3 = @() diag(diag(A));
    D4 = @() I .* dVal;
    t{1}(ii) = timeit(D1,1);
    t{2}(ii) = timeit(D2,1);
    t{3}(ii) = timeit(D3,1);
    t{4}(ii) = timeit(D4,1);
end
maxNumCompThreads('automatic');
%%
figure; hold on
ax = gca;
method{1} = 'sparse(iRow,iCol,dVal,size(A,1),size(A,2))';
method{2} = 'spdiags(B(:,3),0,size(A,1),size(A,2))';
method{3} = 'diag(diag(A))';
method{4} = 'I .* dVal';
plot(matSize.^2,t{1},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{1})
plot(matSize.^2,t{2},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{2})
plot(matSize.^2,t{3},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{3})
plot(matSize.^2,t{4},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{4})
title({'Extracting Sparse Diagonal Matrix from Sparse Matrix';...
    'Generated from Finite Difference of Poisson Equation'})
xlabel('DOF');
ylabel('Time (seconds)')
ax.XScale = 'log';
ax.YScale = 'log';
legend('location','best')