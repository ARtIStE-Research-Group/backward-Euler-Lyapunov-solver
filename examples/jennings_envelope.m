function value = jennings_envelope(t,t1,t2,beta)
%JENNINGS_ENVELOPE Smooth build-up, plateau and exponential decay envelope.

value = zeros(size(t));
idx1 = t < t1;
idx2 = t >= t1 & t <= t2;
idx3 = t > t2;
value(idx1) = (t(idx1)/t1).^2;
value(idx2) = 1;
value(idx3) = exp(-beta*(t(idx3)-t2));
end
