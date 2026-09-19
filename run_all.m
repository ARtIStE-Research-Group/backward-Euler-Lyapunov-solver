%RUN_ALL Execute the automated tests and the principal example.
%
% This script does not require the Control System Toolbox. The package includes
% a complex-Schur Bartels-Stewart Lyapunov solver implemented with base MATLAB.

startup;
run_tests;

fprintf('\nRunning example_01_sdof_modulated ...\n');
example_01_sdof_modulated;
