clear
close all

nDOF = 2.^(4:14);
nHist = 3;

for ii = 1:length(nDOF)
    X = rand(nDOF(ii),nHist);
    F = rand(nDOF(ii),nHist);
    
    XF = (X + F);
    FF = ((F.') * F);
    XFiFF = XF/FF;
    
    P1 = @() (X + F);
    P2 = @() ((F.') * F);
    P3 = @() XF / FF;
    P4 = @() XFiFF * (F.');
    P5 = @() (X + F)/((F.') * F) * (F.');
    
    t{1}(ii) = timeit(P1);
    t{2}(ii) = timeit(P2);
    t{3}(ii) = timeit(P3);
    t{4}(ii) = timeit(P4);
    t{5}(ii) = timeit(P5);
end
%%
figure; hold on
ax = gca;
method{1} = '(X + F)';
method{2} = "((F.') * F)";
method{3} = 'XF / FF';
method{4} = "XFiFF * (F.')";
method{5} = "(X + F)/((F.') * F) * (F.')";
plot(nDOF,t{1},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{1})
plot(nDOF,t{2},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{2})
plot(nDOF,t{3},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{3})
plot(nDOF,t{4},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{4})
plot(nDOF,t{5},'LineWidth',2,'Marker','o','MarkerFaceColor','Auto','DisplayName',method{5})
xlabel('DOF');
ylabel('Time (seconds)')
ax.XScale = 'log';
ax.YScale = 'log';
legend('location','best')