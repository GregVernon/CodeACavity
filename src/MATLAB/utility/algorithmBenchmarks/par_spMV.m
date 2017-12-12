clear
close all
maxNumCompThreads(1);
maxNumCompThreads('automatic');

n = 2:12;
N = 2.^n;
for ii = 1:length(N)
    disp(n(ii))
    A = gallery('poisson',N(ii));
    b = rand(size(A,1),1);
    
    pA = cell(4,1);
    blocks = ceil(linspace(0,size(A,1),5));
    for block = 1:length(blocks)-1
        pA{block} = A(blocks(block)+1:blocks(block+1),:);
    end
    
    tic
    x = cell(length(pA),1);
    parfor block = 1:length(pA)
        x{block} = pA{block}*b;
    end
    t1(ii) = toc;
    x1 = cell2mat(x);
    
    tic
    x2 = A*b;
    t2(ii) = toc;
    
    assert(all(x1==x2))
end
maxNumCompThreads('automatic');

%%
figure; hold on
ax = gca;
method{1} = 'Parallel';
method{2} = 'Serial Multithreaded';

plot(N.^2,t1,'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{1})
plot(N.^2,t2,'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{2})

xlabel('DOF');
ylabel('Time (seconds)')
ax.XScale = 'log';
ax.YScale = 'log';
legend('location','best')

