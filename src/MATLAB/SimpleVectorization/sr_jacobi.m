% function SR_Jacobi(PSI,OMEGA,h,NY,NX,tol)
clear
% Compute relaxation sequence
N = 16;
kmin = sin(pi/2/N)^2;
k = [kmin:kmin:kmin];
k = [k;[k(end) + kmin:kmin:2]'];

% P = 3 scheme
W = [64.66, 6.215, 0.7042];
Q = [1 5 21];

M = sum(Q);
wt = zeros(1,M);
G = ones(size(k));

wt(1) = W(1);
Q(1) = Q(1) - 1;
G = G .* abs(1 - k*wt(1));
counter = 2;
totPerc = 100;
while sum(Q~=0)
    if sum(Q==0) ~= 0
        index = find(Q==0);
        if index == 1
            W = W(2:end);
            Q = Q(2:end);
        elseif index == length(Q)
            W = W(1:end-1);
            Q = Q(1:end-1);
        else
            W(index:end-1) = W(index+1:end);
            W = W(1:end-1);
            Q(index:end-1) = Q(index+1:end);
            Q = Q(1:end-1);
        end
    end
    ww = 1/k(G==max(G));
    dis = abs(W-ww);
    index = find(dis==min(dis));
    wt(counter) = W(index);
    G = G .* abs(1-k*wt(counter));
    Q(index) = Q(index)-1;
    counter = counter + 1;
    perc = sum(Q)/length(wt);
end
% end

%%
clear

% solve(eqn)
% eqn1 = (1-w)*eye(size(J)) + w*J

% solve(eqn1,w)
% inv(D) * (b - R*


