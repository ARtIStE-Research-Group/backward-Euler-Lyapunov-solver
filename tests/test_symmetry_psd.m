function test_symmetry_psd
%TEST_SYMMETRY_PSD Check covariance symmetry and positive semidefiniteness.

A = [0 1; -36 -0.9];
g = [0;1];
B = @(t) 0.10*(1-exp(-0.5*t))^2*(g*g');
opts = lyap_be_options('RelTol',1e-5,'AbsTol',1e-8,'InitialStep',0.01, ...
    'MaxStep',0.25,'CheckPSD',true,'ProjectPSD',false);
sol = lyap_be_adaptive(A,B,[0 20],zeros(2),opts);

symmetryError = 0;
minimumEigenvalue = Inf;
for k = 1:numel(sol.t)
    R = sol.R(:,:,k);
    symmetryError = max(symmetryError,norm(R-R','fro'));
    minimumEigenvalue = min(minimumEigenvalue,min(eig(0.5*(R+R'))));
end
assert(symmetryError < 1e-11,'Symmetry error too large: %.3e',symmetryError);
assert(minimumEigenvalue > -1e-9,'PSD violation too large: %.3e',minimumEigenvalue);
end
