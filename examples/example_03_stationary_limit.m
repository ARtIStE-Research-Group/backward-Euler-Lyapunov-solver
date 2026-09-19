%EXAMPLE_03_STATIONARY_LIMIT Compare the transient result with the algebraic limit.

thisFile = mfilename('fullpath');
packageRoot = fileparts(fileparts(thisFile));
addpath(fullfile(packageRoot,'src'));

A = [0 1; -25 -1.5];
g = [0;1];
B = 0.08*(g*g');
R0 = zeros(2);
opts = lyap_be_options('RelTol',1e-7,'InitialStep',0.01, ...
    'MaxStep',0.25,'CheckPSD',true);
sol = lyap_be_adaptive(A,B,[0 30],R0,opts);

fac = lyap_schur_factor(A);
Rstationary = lyap_schur_solve(fac,B);
Rfinal = sol.R(:,:,end);
relativeError = norm(Rfinal-Rstationary,'fro')/norm(Rstationary,'fro');

fprintf('Relative error with respect to stationary limit: %.3e\n',relativeError);
disp('Stationary covariance:'); disp(Rstationary);
disp('Final transient covariance:'); disp(Rfinal);
