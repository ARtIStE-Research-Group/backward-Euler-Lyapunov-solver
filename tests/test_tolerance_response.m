function test_tolerance_response
%TEST_TOLERANCE_RESPONSE A tighter tolerance should not use fewer accepted steps.

A = [0 1; -9 -0.6];
g = [0;1];
B = @(t) 0.05*(sin(pi*min(t,1)/2).^2)*(g*g');
loose = lyap_be_options('RelTol',1e-3,'InitialStep',0.02, ...
    'MaxStep',0.5,'CheckPSD',false);
tight = lyap_be_options('RelTol',1e-6,'InitialStep',0.02, ...
    'MaxStep',0.5,'CheckPSD',false);
solLoose = lyap_be_adaptive(A,B,[0 8],zeros(2),loose);
solTight = lyap_be_adaptive(A,B,[0 8],zeros(2),tight);
assert(solTight.stats.acceptedSteps >= solLoose.stats.acceptedSteps, ...
    'Tighter tolerance unexpectedly used fewer accepted steps.');
end
