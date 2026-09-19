%EXAMPLE_04_SHEAR_BUILDING Six-storey shear building under modulated base noise.
%
% The input is an illustrative amplitude-modulated white-noise ground
% acceleration. For a colored earthquake model, augment the state with a
% Kanai-Tajimi or Clough-Penzien filter and define the corresponding A and B.

thisFile = mfilename('fullpath');
packageRoot = fileparts(fileparts(thisFile));
addpath(fullfile(packageRoot,'src'));
addpath(fullfile(packageRoot,'examples'));

n = 6;
masses = 2.0e5*ones(n,1);
storyStiffness = linspace(8.0e7,5.0e7,n).';
M = diag(masses);
K = zeros(n);
for story = 1:n
    k = storyStiffness(story);
    K(story,story) = K(story,story)+k;
    if story > 1
        K(story-1,story-1) = K(story-1,story-1)+k;
        K(story,story-1) = K(story,story-1)-k;
        K(story-1,story) = K(story-1,story)-k;
    end
end

% Rayleigh damping: 2% in the first and third modes.
[~,Lambda] = eig(K,M);
omega = sort(sqrt(real(diag(Lambda))));
zetaTarget = 0.02;
rayleigh = [1/(2*omega(1)) omega(1)/2; ...
            1/(2*omega(3)) omega(3)/2] \ [zetaTarget;zetaTarget];
C = rayleigh(1)*M + rayleigh(2)*K;

A = [zeros(n) eye(n); -M\K -M\C];
tau = ones(n,1);
g = [zeros(n,1); -tau];
S0 = 1.0e-2;
phi = @(t) jennings_envelope(t,2.0,14.0,0.25);
B = @(t) 2*pi*S0*phi(t).^2*(g*g');

stateScale = [0.05*ones(n,1); 0.50*ones(n,1)];
opts = lyap_be_options('RelTol',3e-4,'AbsTol',1e-6, ...
    'InitialStep',0.005,'MaxStep',0.20,'StateScale',stateScale, ...
    'CheckPSD',false);
sol = lyap_be_adaptive(A,B,[0 30],zeros(2*n),opts);

D = eye(n);
for floor = 2:n
    D(floor,floor-1) = -1;
end
sigmaDrift = zeros(n,numel(sol.t));
for it = 1:numel(sol.t)
    Ruu = sol.R(1:n,1:n,it);
    Rdrift = D*Ruu*D';
    sigmaDrift(:,it) = sqrt(max(0,diag(Rdrift)));
end

figure('Name','Six-storey shear-building drift response');
plot(sol.t,sigmaDrift.','LineWidth',1.1);
grid on; xlabel('Time'); ylabel('Interstorey-drift standard deviation');
legend(arrayfun(@(i) sprintf('Storey %d',i),1:n,'UniformOutput',false), ...
    'Location','best');
title('Non-stationary interstorey-drift covariance response');

fprintf('Six-storey example: %d accepted steps, %d rejected steps.\n', ...
    sol.stats.acceptedSteps,sol.stats.rejectedSteps);
