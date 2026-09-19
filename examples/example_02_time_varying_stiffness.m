%EXAMPLE_02_TIME_VARYING_STIFFNESS Demonstrate a genuinely time-varying A(t).

thisFile = mfilename('fullpath');
packageRoot = fileparts(fileparts(thisFile));
addpath(fullfile(packageRoot,'src'));
addpath(fullfile(packageRoot,'examples'));

zeta = 0.04;
omega0 = 2*pi*1.2;
A = @(t) [0 1; -(omega0*(1+0.20*sin(0.25*t)))^2, ...
    -2*zeta*omega0*(1+0.20*sin(0.25*t))];
g = [0;1];
S0 = 1.0e-2;
phi = @(t) jennings_envelope(t,1.5,10.0,0.25);
B = @(t) 2*pi*S0*phi(t).^2*(g*g');

opts = lyap_be_options('RelTol',3e-6,'InitialStep',0.01, ...
    'MaxStep',0.25,'StateScale',[0.1;1.0],'CheckPSD',true);
sol = lyap_be_adaptive(A,B,[0 25],zeros(2),opts);

sigmaX = sqrt(max(0,squeeze(sol.R(1,1,:))));
figure('Name','Time-varying state matrix');
plot(sol.t,sigmaX,'LineWidth',1.4); grid on;
xlabel('Time'); ylabel('sigma_x');
title('Response with time-varying stiffness');
