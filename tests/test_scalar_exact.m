function test_scalar_exact
%TEST_SCALAR_EXACT Validate against the analytical scalar solution.

a = -0.8;
b = 2.0;
R0 = 0.3;
tf = 4.0;
opts = lyap_be_options('RelTol',1e-6,'AbsTol',1e-9, ...
    'InitialStep',0.02,'MaxStep',0.2,'CheckPSD',true);
sol = lyap_be_adaptive(a,b,[0 tf],R0,opts);
Rss = -b/(2*a);
Rexact = Rss + (R0-Rss)*exp(2*a*tf);
errorValue = abs(sol.R(:,:,end)-Rexact);
assert(errorValue < 5e-5, ...
    'Scalar exact-solution error is too large: %.3e',errorValue);
end
