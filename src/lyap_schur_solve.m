function [X,relativeResidual] = lyap_schur_solve(fac,Q)
%LYAP_SCHUR_SOLVE Solve A*X + X*A' + Q = 0 from a Schur factorization.
%
%   X = LYAP_SCHUR_SOLVE(FAC,Q) uses FAC returned by LYAP_SCHUR_FACTOR.
%   The algorithm is a complex-Schur Bartels-Stewart method. Columns of the
%   transformed solution are computed in reverse order.

n = fac.n;
if ~(isnumeric(Q) && isequal(size(Q),[n n]) && all(isfinite(Q(:))))
    error('lyap_schur_solve:InvalidQ','Q must be a finite %d-by-%d matrix.',n,n);
end

U = fac.U;
T = fac.T;
Qhat = U' * Q * U;
Y = complex(zeros(n,n));
I = eye(n);

for j = n:-1:1
    rhs = -Qhat(:,j);
    if j < n
        rhs = rhs - Y(:,j+1:n) * conj(T(j,j+1:n)).';
    end
    Mj = T + conj(T(j,j))*I;
    if min(abs(diag(Mj))) < 10*eps*max(1,norm(Mj,1))
        error('lyap_schur_solve:NearlySingular', ...
            ['The Lyapunov operator is singular or nearly singular. ', ...
             'Check whether lambda_i(A)+conj(lambda_j(A)) can vanish.']);
    end
    Y(:,j) = Mj \ rhs;
end

X = U * Y * U';
if isreal(fac.A) && isreal(Q)
    X = real(X);
end
X = 0.5*(X + X');

if nargout > 1
    numerator = norm(fac.A*X + X*fac.A' + Q,'fro');
    denominator = max(1,norm(Q,'fro'));
    relativeResidual = numerator/denominator;
end
end
