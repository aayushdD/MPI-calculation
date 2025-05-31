function O_inf = computeMPI_Zonotope(X, A, W)
    tol=1e-6;
    maxIter=100;
    Omega = X;
    
    for k = 1:maxIter
        OmegaMinusW = pontryDiff(Omega, W);

        PreOmega=inv(A)*OmegaMinusW;
        Omega_new = Omega & PreOmega;

        if isconZonoContainedv2(Omega,Omega_new,tol)
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
