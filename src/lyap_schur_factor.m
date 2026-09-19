function fac = lyap_schur_factor(A)
%LYAP_SCHUR_FACTOR Complex-Schur factorization for Lyapunov solves.
%
%   FAC = LYAP_SCHUR_FACTOR(A) prepares the solution of
%
%       A*X + X*A' + Q = 0
%
%   for one or more right-hand sides Q. The factorization is based on the
%   complex Schur form and requires only base MATLAB.

if ~(isnumeric(A) && ismatrix(A) && size(A,1)==size(A,2))
    error('lyap_schur_factor:SquareMatrix','A must be a square numeric matrix.');
end
if any(~isfinite(A(:)))
    error('lyap_schur_factor:FiniteMatrix','A must contain finite values.');
end

% The complex Schur form avoids the 1x1/2x2 block logic of the real form.
[U,T] = schur(A,'complex');
fac = struct();
fac.A = A;
fac.U = U;
fac.T = T;
fac.n = size(A,1);
fac.eigenvalues = diag(T);
end
