%STARTUP Add the solver, examples and tests to the MATLAB search path.
%
% Run this file from the package root:
%   startup

rootFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(rootFolder, 'src'));
addpath(fullfile(rootFolder, 'examples'));
addpath(fullfile(rootFolder, 'tests'));

fprintf('Adaptive BE Lyapunov Solver paths added.\n');
