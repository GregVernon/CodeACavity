A = [-2 1 0; 1 -2 1; 0 1 -2];
D = eye(size(A)) .* diag(A,0);
R = triu(A,1) + tril(A,-1);
b = zeros(size(A,1),1);

J = inv(D)*b - inv(D)*R;

syms w

W(1) = 2 + sqrt(2);
W(2) = 1;
W(3) = 2 - sqrt(2);
W = sym('w',[3 1]);

for ii = 1:3
    eqn{ii} = (1-W(ii))*eye(size(J)) + W(ii)*J;
end

EQN = sym(1);
for ii = 1:3
    EQN = EQN * eqn{ii};
    % EQN = eqn{1} * eqn{2} * eqn{3}
end
w = solve(EQN);
w = unique(w.w1);
[sw, swidx] = sort(eval(w));
w = w(swidx);