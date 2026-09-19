function M = lyap_eval_matrix(inputValue,t,name)
%LYAP_EVAL_MATRIX Evaluate a constant matrix or a time-dependent callback.

if isa(inputValue,'function_handle')
    M = inputValue(t);
else
    M = inputValue;
end
if ~(isnumeric(M) && ismatrix(M) && all(isfinite(M(:))))
    error('lyap_eval_matrix:InvalidMatrix','%s(t) must return a finite numeric matrix.',name);
end
if ~isfloat(M)
    M = double(M);
end
end
