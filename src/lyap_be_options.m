function opts = lyap_be_options(varargin)
%LYAP_BE_OPTIONS Create or update options for LYAP_BE_ADAPTIVE.
%
%   OPTS = LYAP_BE_OPTIONS() returns the default option structure.
%
%   OPTS = LYAP_BE_OPTIONS('Name',Value,...) overwrites selected defaults.
%
%   OPTS = LYAP_BE_OPTIONS(S) merges the fields of structure S into the
%   defaults. Unknown option names generate an error.
%
% Principal options
%   RelTol              relative tolerance (default 1e-5)
%   AbsTol              absolute tolerance in scaled covariance coordinates
%                       (default 1e-9)
%   InitialStep         first attempted step; [] selects an automatic value
%   MinStep             minimum admissible step; [] selects an automatic value
%   MaxStep             maximum admissible step; [] selects span/10
%   Controller          'P' or 'PI' (default 'PI')
%   StateScale          positive state scale vector s. The covariance is scaled
%                       entrywise by s*s'. Empty means s = ones(m,1).
%   CheckPSD            compute the minimum eigenvalue at accepted steps
%   ProjectPSD          clip negative eigenvalues after acceptance (default false)
%   Verbose             print accepted/rejected-step information
%   OutputFcn           function handle stop = f(t,R,flag)
%
% See also LYAP_BE_ADAPTIVE.

opts = struct();
opts.RelTol             = 1.0e-5;
opts.AbsTol             = 1.0e-9;
opts.InitialStep        = [];
opts.MinStep            = [];
opts.MaxStep            = [];
opts.Safety              = 0.90;
opts.MinFactor           = 0.20;
opts.MaxFactor           = 4.00;
opts.RejectMaxFactor     = 0.80;
opts.Controller          = 'PI';
opts.PIAlpha             = 0.35;  % 0.7/(p+1), p=1
opts.PIBeta              = 0.20;  % 0.4/(p+1), p=1
opts.StateScale          = [];
opts.MaxSteps            = 100000;
opts.SymmetrizeInput     = true;
opts.SymmetrizeSolution  = true;
opts.CheckPSD             = true;
opts.ProjectPSD           = false;
opts.PSDTolerance         = 1.0e-10;
opts.ComputeResidual      = false;
opts.ResidualTolerance    = 1.0e-9;
opts.Verbose              = false;
opts.OutputFcn            = [];

if nargin == 0
    return;
end

if nargin == 1 && isstruct(varargin{1})
    opts = merge_struct(opts, varargin{1});
    opts = validate_options(opts);
    return;
end

if mod(nargin,2) ~= 0
    error('lyap_be_options:NameValuePairs', ...
        'Options must be supplied as name-value pairs or as one structure.');
end

for k = 1:2:nargin
    name = varargin{k};
    value = varargin{k+1};
    if ~(ischar(name) || (isstring(name) && isscalar(name)))
        error('lyap_be_options:InvalidName','Option names must be text.');
    end
    name = char(name);
    names = fieldnames(opts);
    idx = find(strcmpi(name,names),1);
    if isempty(idx)
        error('lyap_be_options:UnknownOption','Unknown option "%s".',name);
    end
    opts.(names{idx}) = value;
end
opts = validate_options(opts);
end

function out = merge_struct(out, in)
names = fieldnames(in);
valid = fieldnames(out);
for k = 1:numel(names)
    idx = find(strcmpi(names{k},valid),1);
    if isempty(idx)
        error('lyap_be_options:UnknownOption','Unknown option "%s".',names{k});
    end
    out.(valid{idx}) = in.(names{k});
end
end

function opts = validate_options(opts)
positiveScalars = {'RelTol','AbsTol','Safety','MinFactor','MaxFactor', ...
    'RejectMaxFactor','PIAlpha','PIBeta','PSDTolerance','ResidualTolerance'};
for k = 1:numel(positiveScalars)
    value = opts.(positiveScalars{k});
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
        error('lyap_be_options:InvalidValue','%s must be a positive finite scalar.',positiveScalars{k});
    end
end
if opts.MinFactor > opts.MaxFactor
    error('lyap_be_options:FactorBounds','MinFactor cannot exceed MaxFactor.');
end
if opts.RejectMaxFactor >= 1
    error('lyap_be_options:RejectFactor','RejectMaxFactor must be smaller than one.');
end
if ~(isnumeric(opts.MaxSteps) && isscalar(opts.MaxSteps) && opts.MaxSteps >= 1)
    error('lyap_be_options:MaxSteps','MaxSteps must be a positive integer.');
end
opts.MaxSteps = floor(opts.MaxSteps);
controller = upper(char(opts.Controller));
if ~ismember(controller,{'P','PI'})
    error('lyap_be_options:Controller','Controller must be ''P'' or ''PI''.');
end
opts.Controller = controller;
logicalFields = {'SymmetrizeInput','SymmetrizeSolution','CheckPSD', ...
    'ProjectPSD','ComputeResidual','Verbose'};
for k = 1:numel(logicalFields)
    if ~(islogical(opts.(logicalFields{k})) || isnumeric(opts.(logicalFields{k})))
        error('lyap_be_options:LogicalOption','%s must be logical.',logicalFields{k});
    end
    opts.(logicalFields{k}) = logical(opts.(logicalFields{k}));
end
if ~isempty(opts.OutputFcn) && ~isa(opts.OutputFcn,'function_handle')
    error('lyap_be_options:OutputFcn','OutputFcn must be empty or a function handle.');
end
end
