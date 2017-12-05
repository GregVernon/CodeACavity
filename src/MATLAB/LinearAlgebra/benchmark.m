clear
n = 3:11;
N = 2.^n;
% N = [20 40 80 160];
nSteps = 1e3;

method = {'Direct', 'Decomposition', 'CG', 'Jacobi'};

f = cell(length(N),length(method));
t = zeros(length(N),length(method));
for ii = 1:length(N)
    for m = 1:length(method)
        f{ii,m} = @() lidCavity(N(ii),1e6,method{m},1e-5,nSteps,false);
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