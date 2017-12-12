clear
close all
ii = 0;
for n = 2:11
    ii = ii+1;
    N(ii) = 2^n;
    
    A1 = gallery('poisson',N(ii));
    B1 = rand(N(ii)^2,1);
    
    B2 = reshape(B1,N(ii),N(ii));
    A2 = rand(size(B2));
    
    gA1 = gpuArray(A1);
    gB1 = gpuArray(B1);
    gA2 = gpuArray(A2);
    gB2 = gpuArray(B2);
    
    f{1} = @() A1*B1;
    f{2} = @() B2(1:end-2,2:end-1) + B2(3:end,2:end-1) + B2(2:end-1,1:end-2) + B2(2:end-1,3:end) - 4*B2(2:end-1,2:end-1);
    f{3} = @() arrayfun(@AX,B2(1:end-2,2:end-1), B2(3:end,2:end-1), B2(2:end-1,1:end-2), B2(2:end-1,3:end), B2(2:end-1,2:end-1));
    f{4} = @() gA1*gB1;
    f{5} = @() gB2(1:end-2,2:end-1) + gB2(3:end,2:end-1) + gB2(2:end-1,1:end-2) + gB2(2:end-1,3:end) - 4*gB2(2:end-1,2:end-1);
    f{6} = @() arrayfun(@AX,gB2(1:end-2,2:end-1), gB2(3:end,2:end-1), gB2(2:end-1,1:end-2), gB2(2:end-1,3:end), gB2(2:end-1,2:end-1));
    
    
    for jj = 1:length(f)
        if ii == 1
            t{jj}(1) = timeit(f{jj});
        else
            if t{jj}(ii-1) <= 1e-3
                try
                    t{jj}(ii) = timeit(f{jj});
                catch
                    t{jj}(ii) = nan;
                end
            else
                t{jj}(ii) = nan;
            end
        end
    end
end
%%
figure; hold on
ax = gca;
plot(N,t{1},'DisplayName','M-V  * CPU');
plot(N,t{2},'DisplayName','Mapped M-V * CPU');
plot(N,t{3},'DisplayName','Mapped F-F arrayfun CPU');
plot(N,t{4},'DisplayName','M-V * GPU');
plot(N,t{5},'DisplayName','Mapped M-V * GPU');
plot(N,t{6},'DisplayName','Mapped F-F arrayfun GPU');
ax.YScale = 'log';
ax.XScale = 'log';
legend