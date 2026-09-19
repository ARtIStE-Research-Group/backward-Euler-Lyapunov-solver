function [R,minEigenvalue,projected] = lyap_make_psd(R,opts)
%LYAP_MAKE_PSD Check and optionally project a covariance matrix onto the PSD cone.

R = 0.5*(R+R');
e = eig(R);
minEigenvalue = min(real(e));
projected = false;
scale = max(1,norm(R,2));
tol = opts.PSDTolerance*scale;

if minEigenvalue < -tol && opts.ProjectPSD
    [V,D] = eig(R);
    d = real(diag(D));
    d(d<0) = 0;
    R = V*diag(d)*V';
    R = real(0.5*(R+R'));
    minEigenvalue = min(real(eig(R)));
    projected = true;
end
end
