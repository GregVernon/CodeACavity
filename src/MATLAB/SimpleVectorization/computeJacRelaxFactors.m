function w = computeJacRelaxFactors(FDM)
% Do not recommend doing for FDM > 4x4
D = eye(size(FDM)) .* diag(FDM,0);
R = triu(FDM,1) + tril(FDM,-1);
b = zeros(size(FDM,1),1);

J = D\b - D\R;

W = sym('w',[size(J,1) 1]);

eqn = cell(length(W),1);
for ii = 1:length(W)
    eqn{ii} = (1-W(ii))*eye(size(J)) + W(ii)*J;
end

EQN = sym(1);
for ii = 1:length(eqn)
    EQN = EQN * eqn{ii};
    % EQN = eqn{1} * eqn{2} * eqn{3}
end
w = solve(EQN);
w = unique(w.w1);
[sw, swidx] = sort(eval(w));
w = w(swidx);
