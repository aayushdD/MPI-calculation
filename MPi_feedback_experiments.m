mu = 4;
ell = 4;
delta = 1;

% K_c: stiffness matrix
Kc = -2*eye(ell) - diag(ones(ell-1,1),1) - diag(ones(ell-1,1),-1);
Kc(1,1) = 1; Kc(end,end) = 1;

% D_c: damping input matrix
Dc = zeros(ell,2);
Dc(1,1) = 1;
Dc(end,2) = -1;

Ac = [zeros(ell), eye(ell);
      -Kc/mu, -delta/mu * eye(ell)];

Bc = [zeros(ell,2);
      Dc/mu];

A = eye(2*ell) + Ac;  % delta = 1
B = Bc;

% Choose LQR weights
% Q = eye(2*ell);
% R = eye(2);

Q = eye(size(A));     % penalize all states equally
R = 0.1 * eye(size(B,2));  % small input penalty
K = -lqr(A, B, Q, R);

X=zono(eye(ell*2),zeros(ell*2,1));
U=zono(eye(2),zeros(2,1));
% W=zono(zeros(ell*2),zeros(ell*2,1));

MPI=computeMPI_Zonotope_feedback(Zi, U, A, B, [], K);