%EXAMPLE_01_SDOF_MODULATED Linear SDOF under amplitude-modulated white noise.
%
% State: z = [relative displacement; relative velocity].
% The example illustrates covariance histories and adaptive step selection.

thisFile = mfilename('fullpath');
packageRoot = fileparts(fileparts(thisFile));
addpath(fullfile(packageRoot,'src'));
addpath(fullfile(packageRoot,'examples'));

mass = 1.0;
frequencyHz = 1.5;
omega = 2*pi*frequencyHz;
zeta = 0.05;
A = [0 1; -omega^2 -2*zeta*omega];
g = [0; 1/mass];
S0 = 2.0e-2;
phi = @(t) jennings_envelope(t,2.0,12.0,0.30);
B = @(t) 2*pi*S0*phi(t).^2*(g*g');

opts = lyap_be_options( ...
    'RelTol',1e-4, ...
    'AbsTol',1e-6, ...
    'InitialStep',0.01, ...
    'MaxStep',0.50, ...
    'StateScale',[0.10; 1.00], ...
    'Controller','PI', ...
    'CheckPSD',true, ...
    'Verbose',false);

sol = lyap_be_adaptive(A,B,[0 30],zeros(2),opts);

sigmaX = sqrt(max(0,squeeze(sol.R(1,1,:))));
sigmaV = sqrt(max(0,squeeze(sol.R(2,2,:))));

figure('Name','Adaptive Lyapunov response');
plot(sol.t,sigmaX,'LineWidth',1.4); hold on;
plot(sol.t,sigmaV,'LineWidth',1.4);
grid on; xlabel('Time'); ylabel('Standard deviation');
legend('sigma_x','sigma_v','Location','best');
title('Non-stationary covariance response');

figure('Name','Accepted and rejected steps');
accepted = sol.attempts.accepted;
semilogy(sol.attempts.t(accepted),sol.attempts.h(accepted),'o-','LineWidth',1.0); hold on;
if any(~accepted)
    semilogy(sol.attempts.t(~accepted),sol.attempts.h(~accepted),'x','MarkerSize',8,'LineWidth',1.2);
    legend('Accepted','Rejected','Location','best');
else
    legend('Accepted','Location','best');
end
grid on; xlabel('Step starting time'); ylabel('Attempted step size');
title('Adaptive step-size history');

fprintf('\nExample 1 completed.\n');
fprintf('Accepted steps : %d\n',sol.stats.acceptedSteps);
fprintf('Rejected steps : %d\n',sol.stats.rejectedSteps);
fprintf('Lyapunov solves: %d\n',sol.stats.lyapunovSolves);
fprintf('Schur factors  : %d\n',sol.stats.schurFactorizations);
