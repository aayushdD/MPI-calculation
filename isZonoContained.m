function [isContained, diagnostics] = isZonoContained(Y, X)
% Check if constrained zonotope Y is contained in X using strong duality (Equation 12)
% via a linear program using Gurobi

% Extract zonotope data
cX = X.c; GX = X.G; AX = X.A; bX = X.b;
cY = Y.c; GY = Y.G; AY = Y.A; bY = Y.b;

n = length(cX);
NX = size(GX, 2);      % dim of xi_X
MX = size(AX, 1);      % num of constraints in X
NY = size(GY, 2);      % dim of eta_Y
MY = size(AY, 1);      % num of constraints in Y

% Fix AY dimension if transposed
if size(AY, 2) ~= MY
    AY = AY';  % Ensure AY is MY x nBeta
end

% Variable counts
nAlpha = n;
nBeta = MY;
nXi = NX;
nZ = NY + MY;  % for 1-norm: |GY'*alpha + AY'*beta| <= 1

nVars = nAlpha + nBeta + nXi + nZ;

% Variable bounds
lb = [-inf(nAlpha + nBeta + nXi, 1); zeros(nZ,1)];
ub = inf(nVars,1);

% Objective: minimize 1 + alpha' * (cY - GX * xi - cX) - beta' * bY
obj = zeros(nVars, 1);
obj(1:nAlpha) = (cY - cX);                 % alpha term
obj(nAlpha+1:nAlpha+nBeta) = -bY;          % beta term

% Equality: AX * xi == bX
Aeq1 = zeros(MX, nVars);
Aeq1(:, nAlpha + nBeta + 1 : nAlpha + nBeta + nXi) = AX;
beq1 = bX;

% 1-norm constraint: ||GY'*alpha + AY'*beta||_1 ≤ 1
Aineq1 = zeros(2*nZ, nVars);
bineq1 = zeros(2*nZ, 1);

% Define u = GY' * alpha + AY' * beta and encode with auxiliary z variables
for i = 1:nZ
    if i <= NY
        % z_i >= GY(:,i)' * alpha
        Aineq1(2*i-1, 1:nAlpha) = -GY(:,i)';
        Aineq1(2*i-1, nAlpha + nBeta + nXi + i) = -1;

        % z_i >= -GY(:,i)' * alpha
        Aineq1(2*i, 1:nAlpha) = GY(:,i)';
        Aineq1(2*i, nAlpha + nBeta + nXi + i) = -1;
    else
        j = i - NY;
        % z_i >= AY(j,:)*beta
        Aineq1(2*i-1, nAlpha+1:nAlpha+nBeta) = -AY(j,:);
        Aineq1(2*i-1, nAlpha + nBeta + nXi + i) = -1;

        % z_i >= -AY(j,:)*beta
        Aineq1(2*i, nAlpha+1:nAlpha+nBeta) = AY(j,:);
        Aineq1(2*i, nAlpha + nBeta + nXi + i) = -1;
    end
end

% Sum(z) ≤ 1
Aineq2 = zeros(1, nVars);
Aineq2(1, nAlpha + nBeta + nXi + 1 : end) = 1;
bineq2 = 1;

% Combine all constraints
model.A = sparse([Aeq1; Aineq1; Aineq2]);
model.rhs = [beq1; bineq1; bineq2];
model.sense = [repmat('=', size(Aeq1,1),1); repmat('<', size(Aineq1,1),1); '<'];
model.modelsense = 'min';
model.obj = obj;
model.lb = lb;
model.ub = ub;

params.outputflag = 0;
result = gurobi(model, params);

% Output
if isfield(result, 'objval') && strcmp(result.status, 'OPTIMAL')
    isContained = result.objval >= -0.5;
    diagnostics.optval = result.objval;
else
    isContained = false;
    diagnostics.optval = NaN;
end
