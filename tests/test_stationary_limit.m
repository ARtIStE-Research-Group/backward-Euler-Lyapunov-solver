function test_stationary_limit
%TEST_STATIONARY_LIMIT Check convergence to the algebraic Lyapunov solution.

A = [0 1; -16 -1.2];
g = [0;1];
B = 0.05*(g*g');
opts = lyap_be_options('RelTol',1e-5,'AbsTol',1e-8,'InitialStep',0.01, ...
    'MaxStep',0.2,'CheckPSD',true);
sol = lyap_be_adaptive(A,B,[0 35],zeros(2),opts);
Rinf = lyap_schur_solve(lyap_schur_factor(A),B);
rel = norm(sol.R(:,:,end)-Rinf,'fro')/norm(Rinf,'fro');
assert(rel < 2e-4,'Stationary-limit relative error is too large: %.3e',rel);
end
