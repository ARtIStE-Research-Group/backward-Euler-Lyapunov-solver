function err = lyap_be_error_norm(E,Rhigh,Rold,opts)
%LYAP_BE_ERROR_NORM Weighted RMS Frobenius norm for step-doubling error.
%
% The optional StateScale vector s transforms covariance entries according to
%   R_scaled(i,j) = R(i,j)/(s(i)*s(j)).
% AbsTol and RelTol are then applied componentwise in scaled coordinates.

m = size(E,1);
if isempty(opts.StateScale)
    s = ones(m,1);
else
    s = opts.StateScale(:);
    if numel(s) ~= m || any(~isfinite(s)) || any(s <= 0)
        error('lyap_be_error_norm:StateScale', ...
            'StateScale must contain %d positive finite entries.',m);
    end
end
S = s*s.';
Es = E./S;
Rhs = Rhigh./S;
Ros = Rold./S;
weights = opts.AbsTol + opts.RelTol.*max(abs(Rhs),abs(Ros));
ratio = Es./weights;
err = norm(ratio,'fro')/sqrt(numel(ratio));
end
