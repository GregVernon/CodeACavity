clear
close all
clc
figure; hold on; 
N = [10:1:100];
%%
cd Native
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    %     N = 2^n;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(1,kk) = timeit(FCN);
    if kk == 1
        p{1} = plot(N(1:kk),t(1,1:kk),'DisplayName','Native');
        legend(p{1});
    else
        p{1}.XData = N(1:kk);
        p{1}.YData = t(1,1:kk);
    end
    drawnow
end
cd ..
%%
cd Built-in
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(2,kk) = timeit(FCN);
    if kk == 1
        p{2} = plot(N(1:kk),t(2,1:kk),'DisplayName','Built-in');
    else
        p{2}.XData = N(1:kk);
        p{2}.YData = t(2,1:kk);
    end
    drawnow
end
cd ..

%%
cd SimpleVectorization
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(3,kk) = timeit(FCN);
    if kk == 1
        p{3} = plot(N(1:kk),t(3,1:kk),'DisplayName','SimpleVectorization');
    else
        p{3}.XData = N(1:kk);
        p{3}.YData = t(3,1:kk);
    end
    drawnow
end
cd ..
%%
% cd SimpleVectorization_GPU
% FCN = @() lidCavity(N(ii),10,100,1e-5,false);
% t(4) = timeit(FCN);
% cd ..
%%
cd AdvancedVectorization_GPU
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(5,kk) = timeit(FCN);
    if kk == 1
        p{5} = plot(N(1:kk),t(5,1:kk),'DisplayName','AdvancedVectorization GPU');
    else
        p{5}.XData = N(1:kk);
        p{5}.YData = t(5,1:kk);
    end
    drawnow
end
cd ..
%%
cd LinearAlgebra
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(6,kk) = timeit(FCN);
    if kk == 1
        p{6} = plot(N(1:kk),t(6,1:kk),'DisplayName','LinearAlgebra');
    else
        p{6}.XData = N(1:kk);
        p{6}.YData = t(6,1:kk);
    end
    drawnow
end
cd ..
%%
cd LinearAlgebra_GPU
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(7,kk) = timeit(FCN);
    if kk == 1
        p{7} = plot(N(1:kk),t(7,1:kk),'DisplayName','LinearAlgebra GPU');
    else
        p{7}.XData = N(1:kk);
        p{7}.YData = t(7,1:kk);
    end
    drawnow
end
cd ..
%%
cd LinearAlgebra_CPU_GPU
kk = 0;
for ii = 1:length(N)
    kk = kk+1;
    FCN = @() lidCavity(N(ii),10,100,1e-5,false);
    t(8,kk) = timeit(FCN);
    if kk == 1
        p{8} = plot(N(1:kk),t(8,1:kk),'DisplayName','LinearAlgebra CPU GPU');
    else
        p{8}.XData = N(1:kk);
        p{8}.YData = t(8,1:kk);
    end
    drawnow
end
cd ..
