function sol = lyap_be_adaptive(Ainput,Binput,tspan,R0,options)
%LYAP_BE_ADAPTIVE Variable-step backward-Euler integration of a differential Lyapunov equation.
%
%   SOL = LYAP_BE_ADAPTIVE(A,B,TSPAN,R0,OPTS) integrates
%
%       dR/dt = A(t)R + R A(t)' + B(t),    R(t0)=R0,
%
%   on TSPAN = [t0 tf]. A and B may be constant matrices or function handles
%   A(t), B(t). The method uses backward Euler and step doubling. Every
%   implicit substep is recast as an algebraic Lyapunov equation and solved by
%   a reusable complex-Schur Bartels-Stewart factorization.
%
% Output fields
%   sol.t                  accepted output times, 1-by-N
%   sol.R                  covariance history, m-by-m-by-N
%   sol.h                  accepted step sizes
%   sol.err                accepted normalized local-error estimates
%   sol.minEigenvalue      minimum eigenvalue at each output time (if enabled)
%   sol.attempts           attempted-step history, including rejections
%   sol.stats              counters and termination information
%   sol.options            actual option structure
%
% Example
%   A = [0 1; -100 -1];
%   g = [0;1]; S0 = 0.01;
%   phi = @(t) min(1,t/2).*exp(-0.08*max(t-12,0));
%   B = @(t) 2*pi*S0*phi(t)^2*(g*g');
%   opts = lyap_be_options('RelTol',1e-6,'InitialStep',0.01);
%   sol = lyap_be_adaptive(A,B,[0 30],zeros(2),opts);
%
% See also LYAP_BE_OPTIONS, LYAP_SCHUR_FACTOR, LYAP_SCHUR_SOLVE.

if nargin < 5 || isempty(options)
    opts = lyap_be_options();
elseif isstruct(options)
    opts = lyap_be_options(options);
else
    error('lyap_be_adaptive:Options','The fifth argument must be an option structure.');
end

if ~(isnumeric(tspan) && numel(tspan)==2 && all(isfinite(tspan)))
    error('lyap_be_adaptive:Tspan','tspan must be [t0 tf].');
end
t0 = tspan(1); tf = tspan(2);
if tf <= t0
    error('lyap_be_adaptive:TspanOrder','tf must be greater than t0.');
end

A0 = lyap_eval_matrix(Ainput,t0,'A');
B0 = lyap_eval_matrix(Binput,t0,'B');
m = size(A0,1);
if size(A0,2) ~= m
    error('lyap_be_adaptive:A','A must be square.');
end
if ~isequal(size(B0),[m m])
    error('lyap_be_adaptive:B','B must have the same dimensions as A.');
end
if ~(isnumeric(R0) && isequal(size(R0),[m m]) && all(isfinite(R0(:))))
    error('lyap_be_adaptive:R0','R0 must be a finite matrix with the same dimensions as A.');
end
if ~isfloat(R0)
    R0 = double(R0);
end
if opts.SymmetrizeInput
    B0 = 0.5*(B0+B0');
    R0 = 0.5*(R0+R0');
end
if isempty(opts.StateScale)
    opts.StateScale = ones(m,1);
else
    s = opts.StateScale(:);
    if numel(s) ~= m || any(~isfinite(s)) || any(s<=0)
        error('lyap_be_adaptive:StateScale','StateScale must have %d positive entries.',m);
    end
    opts.StateScale = s;
end

span = tf-t0;
if isempty(opts.MaxStep)
    opts.MaxStep = span/10;
end
if isempty(opts.MinStep)
    opts.MinStep = max(100*eps(max(1,max(abs(tspan)))),1e-12*span);
end
if isempty(opts.InitialStep)
    dynamicStep = 0.1/max(1,norm(A0,1));
    opts.InitialStep = min([span/100,dynamicStep,opts.MaxStep]);
end
if ~(opts.MinStep>0 && opts.MaxStep>=opts.MinStep && opts.InitialStep>0)
    error('lyap_be_adaptive:StepBounds','Invalid InitialStep, MinStep or MaxStep.');
end

h = min(opts.MaxStep,max(opts.MinStep,opts.InitialStep));
t = t0;
R = R0;
constantA = isnumeric(Ainput);
constantB = isnumeric(Binput);
if constantA
    Aconst = A0;
end
if constantB
    Bconst = B0;
end

capacity = 1024;
tout = zeros(1,capacity);
Rout = zeros(m,m,capacity,class(R0));
hout = zeros(1,capacity);
errout = zeros(1,capacity);
minEigOut = nan(1,capacity);
tout(1)=t; Rout(:,:,1)=R;
if opts.CheckPSD
    minEigOut(1)=min(real(eig(0.5*(R+R'))));
end
nAccepted = 0;
nRejected = 0;
nAttempts = 0;
nFactorizations = 0;
nLyapSolves = 0;
nProjected = 0;
maxResidual = 0;
previousError = 1;
warnedPSD = false;
terminatedByOutputFcn = false;

attemptTime = zeros(1,capacity);
attemptH = zeros(1,capacity);
attemptErr = zeros(1,capacity);
attemptAccepted = false(1,capacity);

cache.h = NaN;
cache.full = [];
cache.half = [];

if ~isempty(opts.OutputFcn)
    stop = opts.OutputFcn(t,R,'init');
    if isempty(stop), stop = false; end
    if stop
        terminatedByOutputFcn = true;
    end
end

while t < tf && ~terminatedByOutputFcn
    if nAttempts >= opts.MaxSteps
        error('lyap_be_adaptive:MaxSteps','Maximum number of attempted steps exceeded.');
    end
    h = min(h,tf-t);
    nAttempts = nAttempts+1;

    tmid = t+0.5*h;
    tend = t+h;
    if constantA
        Amid = Aconst; Aend = Aconst;
    else
        Amid = lyap_eval_matrix(Ainput,tmid,'A');
        Aend = lyap_eval_matrix(Ainput,tend,'A');
    end
    if ~isequal(size(Amid),[m m]) || ~isequal(size(Aend),[m m])
        error('lyap_be_adaptive:ChangingDimension','A(t) cannot change dimensions.');
    end
    if constantB
        Bmid = Bconst; Bend = Bconst;
    else
        Bmid = lyap_eval_matrix(Binput,tmid,'B');
        Bend = lyap_eval_matrix(Binput,tend,'B');
    end
    if ~isequal(size(Bmid),[m m]) || ~isequal(size(Bend),[m m])
        error('lyap_be_adaptive:ChangingDimension','B(t) cannot change dimensions.');
    end
    if opts.SymmetrizeInput
        Bmid=0.5*(Bmid+Bmid'); Bend=0.5*(Bend+Bend');
    end

    sameCachedStep = constantA && isfinite(cache.h) && ...
        abs(h-cache.h) <= 10*eps(max(1,abs(h)));
    if sameCachedStep
        facFull = cache.full;
        facHalf1 = cache.half;
        facHalf2 = cache.half;
    elseif constantA
        facFull = lyap_schur_factor(Aconst-eye(m)/(2*h));
        facHalf1 = lyap_schur_factor(Aconst-eye(m)/h);
        facHalf2 = facHalf1;
        cache.h = h; cache.full = facFull; cache.half = facHalf1;
        nFactorizations = nFactorizations+2;
    else
        facFull = lyap_schur_factor(Aend-eye(m)/(2*h));
        facHalf1 = lyap_schur_factor(Amid-eye(m)/h);
        facHalf2 = lyap_schur_factor(Aend-eye(m)/h);
        nFactorizations = nFactorizations+3;
    end

    Qfull = Bend + R/h;
    if opts.ComputeResidual
        [Rfull,res1] = lyap_schur_solve(facFull,Qfull);
    else
        Rfull = lyap_schur_solve(facFull,Qfull); res1=0;
    end
    Qhalf1 = Bmid + 2*R/h;
    if opts.ComputeResidual
        [Rhalf1,res2] = lyap_schur_solve(facHalf1,Qhalf1);
    else
        Rhalf1 = lyap_schur_solve(facHalf1,Qhalf1); res2=0;
    end
    Qhalf2 = Bend + 2*Rhalf1/h;
    if opts.ComputeResidual
        [Rhalf2,res3] = lyap_schur_solve(facHalf2,Qhalf2);
    else
        Rhalf2 = lyap_schur_solve(facHalf2,Qhalf2); res3=0;
    end
    nLyapSolves = nLyapSolves+3;
    maxResidual = max([maxResidual,res1,res2,res3]);
    if opts.ComputeResidual && max([res1,res2,res3]) > opts.ResidualTolerance
        warning('lyap_be_adaptive:Residual', ...
            'A Lyapunov solve exceeded the requested residual tolerance.');
    end

    if opts.SymmetrizeSolution
        Rfull=0.5*(Rfull+Rfull');
        Rhalf2=0.5*(Rhalf2+Rhalf2');
    end
    E = Rhalf2-Rfull;
    err = lyap_be_error_norm(E,Rhalf2,R,opts);
    if ~isfinite(err)
        err = Inf;
    end

    if nAttempts > numel(attemptTime)
        attemptTime = [attemptTime zeros(1,capacity)]; %#ok<AGROW>
        attemptH = [attemptH zeros(1,capacity)]; %#ok<AGROW>
        attemptErr = [attemptErr zeros(1,capacity)]; %#ok<AGROW>
        attemptAccepted = [attemptAccepted false(1,capacity)]; %#ok<AGROW>
    end
    attemptTime(nAttempts)=t;
    attemptH(nAttempts)=h;
    attemptErr(nAttempts)=err;

    if err <= 1
        nAccepted=nAccepted+1;
        attemptAccepted(nAttempts)=true;
        t=tend;
        R=Rhalf2;
        minEig=NaN;
        if opts.CheckPSD
            [R,minEig,projected]=lyap_make_psd(R,opts);
            if projected
                nProjected=nProjected+1;
            end
            scale=max(1,norm(R,2));
            if minEig < -opts.PSDTolerance*scale && ~opts.ProjectPSD && ~warnedPSD
                warning('lyap_be_adaptive:PSD', ...
                    ['A significant negative covariance eigenvalue was detected. ', ...
                     'Inspect tolerances, input B(t), and model stability.']);
                warnedPSD=true;
            end
        end
        outIndex=nAccepted+1;
        if outIndex > numel(tout)
            tout=[tout zeros(1,capacity)]; %#ok<AGROW>
            Rout=cat(3,Rout,zeros(m,m,capacity,class(R0))); %#ok<AGROW>
            hout=[hout zeros(1,capacity)]; %#ok<AGROW>
            errout=[errout zeros(1,capacity)]; %#ok<AGROW>
            minEigOut=[minEigOut nan(1,capacity)]; %#ok<AGROW>
        end
        tout(outIndex)=t; Rout(:,:,outIndex)=R;
        hout(nAccepted)=h; errout(nAccepted)=err;
        minEigOut(outIndex)=minEig;

        if opts.Verbose
            fprintf('ACCEPT t=%12.5g  h=%10.3e  err=%10.3e\n',t,h,err);
        end
        if ~isempty(opts.OutputFcn)
            stop = opts.OutputFcn(t,R,'step');
            if isempty(stop), stop = false; end
            if stop
                terminatedByOutputFcn=true;
            end
        end

        if err == 0
            factor=opts.MaxFactor;
        elseif strcmp(opts.Controller,'PI') && nAccepted>1
            factor=opts.Safety*err^(-opts.PIAlpha)*previousError^(opts.PIBeta);
        else
            factor=opts.Safety*err^(-0.5);
        end
        factor=min(opts.MaxFactor,max(opts.MinFactor,factor));
        previousError=max(err,realmin);
        h=min(opts.MaxStep,max(opts.MinStep,h*factor));
    else
        nRejected=nRejected+1;
        if opts.Verbose
            fprintf('REJECT t=%12.5g  h=%10.3e  err=%10.3e\n',t,h,err);
        end
        if isfinite(err) && err>0
            factor=opts.Safety*err^(-0.5);
        else
            factor=opts.MinFactor;
        end
        factor=min(opts.RejectMaxFactor,max(opts.MinFactor,factor));
        hnew=max(opts.MinStep,h*factor);
        if h <= opts.MinStep*(1+10*eps) && err>1
            error('lyap_be_adaptive:MinimumStep', ...
                'The error tolerance cannot be met at MinStep.');
        end
        h=hnew;
    end
end

if ~isempty(opts.OutputFcn)
    opts.OutputFcn(t,R,'done');
end

N=nAccepted+1;
sol=struct();
sol.t=tout(1:N);
sol.R=Rout(:,:,1:N);
sol.h=hout(1:nAccepted);
sol.err=errout(1:nAccepted);
if opts.CheckPSD
    sol.minEigenvalue=minEigOut(1:N);
else
    sol.minEigenvalue=[];
end
sol.attempts=struct();
sol.attempts.t=attemptTime(1:nAttempts);
sol.attempts.h=attemptH(1:nAttempts);
sol.attempts.err=attemptErr(1:nAttempts);
sol.attempts.accepted=attemptAccepted(1:nAttempts);
sol.stats=struct();
sol.stats.acceptedSteps=nAccepted;
sol.stats.rejectedSteps=nRejected;
sol.stats.attemptedSteps=nAttempts;
sol.stats.lyapunovSolves=nLyapSolves;
sol.stats.schurFactorizations=nFactorizations;
sol.stats.psdProjections=nProjected;
sol.stats.maxRelativeLyapunovResidual=maxResidual;
sol.stats.terminatedByOutputFcn=terminatedByOutputFcn;
sol.stats.finalTime=t;
sol.options=opts;
end
