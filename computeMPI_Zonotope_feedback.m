function O_inf = computeMPI_Zonotope_feedback(X, U, A, B, W, K)
    % Inputs:
    %   X, U, W - zonotopes with fields: c, G
    %   A, B    - system matrices
    %   W       - disturbance zonotope
    %   K       - feedback gain (optional)
    %
    % If B or K are empty, assumes autonomous or open-loop system.

    tol = 1e-6;
    maxIter = 100;

    % Decide system type
    hasInput = ~isempty(B);
    useFeedback = ~isempty(K);
    hasdisturbance=~isempty(W);

    if ~hasInput
        % Autonomous system
        A_cl = A;
        Omega = X;
        fprintf('Autonomous system: using Omega_0 = X\n');
    elseif useFeedback
        % Closed-loop feedback system
        A_cl = A + B * K;

        % Extract matrices
        c_X = X.c;
        G_X = X.G;
        c_U = U.c;
        G_U = U.G;

        % Kx ∈ U ⇒ K G_X * xi = c_U - K c_X
        F_bar = [K * G_X , -G_U];
        theta_bar = c_U - K * c_X;

        % Build constrained zonotope
        X_bar = conZono([G_X, zeros(size(G_X,1),size(G_U,2))],c_X, F_bar, theta_bar);
        Omega = X_bar;
        fprintf('Feedback system: using Omega_0 = X ∩ {x : Kx ∈ U}\n');
    else
        % Open-loop with input but no feedback
        A_cl = A;
        Omega = X;
        fprintf('Open-loop system with input: using Omega_0 = X\n');
    end

    for k = 1:maxIter
        if ~hasdisturbance
            PreOmega=inv(A_cl)*Omega;
        else
            OmegaMinusW = pontryDiff(Omega, W);
            PreOmega=inv(A_cl)*OmegaMinusW;
        end

        Omega_new = Omega & PreOmega;

        if isconZonoContainedv2(Omega, Omega_new, tol)
            Omega = Omega_new;
            fprintf('Algorithm converged after %d iterations.\n', k);
            break;
        end
        Omega = Omega_new;
    end

    O_inf = Omega;

    if k == maxIter
        warning('Maximum iterations reached before convergence.');
    end
end