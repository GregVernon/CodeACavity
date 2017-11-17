clear
close all
n = 3:6;
N = 2.^n;
% N = [20 40 80 160];
nSteps = 1e2;

method = {'Jacobi', 'Weighted-Jacobi'};

f = cell(length(N),length(method));
t = zeros(length(N),length(method));
for ii = 1:length(N)
    for m = 1:length(method)
        f{ii,m} = @() lidCavity(N(ii),inf,method{m},1e-5,nSteps,false);
        t(ii,m) = timeit(f{ii,m});
        disp(num2str(N(ii))+"x"+num2str(N(ii)) + " via " + method{m} + ": " + num2str(t(ii,m)) + " seconds to compute " + num2str(nSteps) + " timesteps")
        pause(0.1)
    end
end
%%
figure;
hold on
for m = 1:length(method)
    plot((N-2).^2, t(:,m),'DisplayName',method{m})
end
legend('Location','Northwest')
ax1 = gca;
ax1.XScale = 'log';
ax1.YScale = 'log';
xlabel('# DOF')
ylabel('Time (s)')
title("Solving " + num2str(nSteps) + " Timesteps of Lid-Cavity Flow")