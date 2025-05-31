function isContained = isconZonoContainedv2(Z1, Z2, tolerance)

    if nargin < 3
        tolerance = 1e-6;  % default small buffer if not provided
    end

    c1 = Z1.c; G1 = Z1.G; F1 = Z1.A; theta1 = Z1.b;
    c2 = Z2.c; G2 = Z2.G; F2 = Z2.A; theta2 = Z2.b;

    [~, D1] = size(G1);
    [~, D2] = size(G2);
    q1 = size(F1, 1);
    q2 = size(F2, 1);

    %% Optimization variables
    Gamma = sdpvar(D2, D1, 'full');
    beta = sdpvar(D2, 1);
    Pi = sdpvar(q2, q1);

    Constraints = [];

    % (19a) with tolerance buffer
    Constraints = [Constraints, norm(c2 - c1 - G2 * beta, Inf) <= tolerance];
    Constraints = [Constraints, norm(G1 - G2 * Gamma, Inf) <= tolerance];

    % (19b) equality constraints remain strict
    Constraints = [Constraints, Pi * F1 == F2 * Gamma];
    Constraints = [Constraints, Pi * theta1 == theta2 + F2 * beta];

    % (19c) bounding the image of zonotope
    Constraints = [Constraints, sum(abs(Gamma), 2) + abs(beta) <= 1 + tolerance];
    Constraints = [Constraints, Pi >= 0];

    %% Solve feasibility problem
    options = sdpsettings('solver', 'gurobi', 'verbose', 0);
    diagnostics = optimize(Constraints, [], options);

    isContained = diagnostics.problem == 0;
end
