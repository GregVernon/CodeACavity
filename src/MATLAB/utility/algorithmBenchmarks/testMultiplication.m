clear
close all
ii = 0;
for n = 2:2:22
    ii = ii+1;
    N(ii) = 2^n;
    
    A1 = gallery('poisson',sqrt(N(ii)));
    B1 = rand(N(ii),1);
    
    B2 = reshape(B1,sqrt(N(ii)),sqrt(N(ii)));
    A2 = rand(size(B2));
    
    gA1 = gpuArray(A1);
    gB1 = gpuArray(B1);
    gA2 = gpuArray(A2);
    gB2 = gpuArray(B2);
    
    f{1} = @() A1*B1;
    f{2} = @() A2.*B2;
    f{3} = @() arrayfun(@times,A2,B2);
    f{4} = @() bsxfun(@times,A2,B2);
    f{5} = @() gA1*gB1;
    f{6} = @() gA2.*gB2;
    f{7} = @() arrayfun(@times,gA2,gB2);
    f{8} = @() bsxfun(@times,gA2,gB2);
    
    for jj = 1:length(f)
        if ii == 1
            t{jj}(1) = timeit(f{jj});
        else
            if t{jj}(ii-1) <= 1
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
plot(N,t{2},'DisplayName','F-F .* CPU');
plot(N,t{3},'DisplayName','F-F arrayfun CPU');
plot(N,t{4},'DisplayName','F-F bsxfun CPU');
plot(N,t{5},'DisplayName','M-V * GPU');
plot(N,t{6},'DisplayName','F-F .* GPU');
plot(N,t{7},'DisplayName','F-F arrayfun GPU');
plot(N,t{8},'DisplayName','F-F bsxfun GPU');
ax.YScale = 'log';
ax.XScale = 'log';
legend