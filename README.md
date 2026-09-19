# Backward-Euler Matrix Lyapunov eq. Solver

This package integrates the non-stationary differential Lyapunov equation

    dR/dt = A(t) R + R A(t)' + B(t)

using backward Euler, step doubling, and proportional/PI step-size control.
Each implicit substep is transformed into an algebraic Lyapunov equation.
A complex-Schur Bartels-Stewart solver is included, so the package does not
require MATLAB's Control System Toolbox.

## Quick start

```matlab
cd backward-Euler-Lyapunov-solver
startup
run_tests
example_01_sdof_modulated
```

## Minimal use

```matlab
A = [0 1; -100 -1];
g = [0;1];
S0 = 0.01;
phi = @(t) min(1,t/2).*exp(-0.08*max(t-12,0));
B = @(t) 2*pi*S0*phi(t)^2*(g*g');

opts = lyap_be_options('RelTol',1e-4, ...
                       'InitialStep',0.01, ...
                       'MaxStep',0.5, ...
                       'StateScale',[0.1;1.0]);
sol = lyap_be_adaptive(A,B,[0 30],zeros(2),opts);

sigmaX = sqrt(squeeze(sol.R(1,1,:)));
plot(sol.t,sigmaX)
```

## Package contents

- `src/`: numerical solver and Schur-based Lyapunov routines.
- `examples/`: constant- and time-varying-state-matrix examples.
- `tests/`: regression and validation tests.
- `manual/Backward_Euler_Matrix_Lyapunov_eq__Solver_English_Manual.pdf`: theory and user manual (English).

## License

Released under the MIT License; see `LICENSE`.

## Citation

If you use this software in research or professional work, please cite the
corresponding software release.

Suggested citation format:

> Marano, G. C., Datta, G., Sardone, L. (2026). *Backward-Euler Matrix
> Lyapunov eq. Solver*. Software repository, ARTISTE Research Group,
> Politecnico di Torino. DOI: 10.5281/zenodo.22842946

See also `CITATION.cff` for the machine-readable citation metadata.

## Important qualification

The implementation is a research software prototype. The supplied tests are
necessary but not sufficient for safety-critical engineering use. Validate the
model, tolerances, units, and output against independent calculations.
