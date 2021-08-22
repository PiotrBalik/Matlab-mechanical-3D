function [x,history,fval,exitFlag,output,population,scores] = gamultiobj2(fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options)
%GAMULTIOBJ Multiobjective optimization using genetic algorithm.
%   GAMULTIOBJ attempts to solve multiobjective problems of the form:
%       min F(X)  subject to:  A*X <= b, Aeq*X = beq (linear constraints)
%        X                     lb <= X <= ub (bound constraints)
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS) finds a local Pareto set X of the
%   objective functions defined in FITNESSFCN. NVARS is the dimension of
%   the optimization problem (number of decision variables). FITNESSFCN
%   accepts a vector X of size 1-by-NVARS and returns a vector of
%   size 1-by-numberOfObjectives evaluated at a decision variable. X is
%   a matrix with NVARS columns. The number of rows in X is the same as the
%   number of Pareto solutions. All solutions in a Pareto set are equally
%   optimal, and it is up to the designer to select a solution in the Pareto
%   set depending on the application.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b) finds a local Pareto set X of the
%   objective functions defined in FITNESSFCN, subject to the linear
%   inequalities A*X <= B. Linear constraints are supported only for
%   default PopulationType option ('doubleVector'); other population types
%   e.g., 'bitString' and 'custom' are not supported.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq) finds a local Pareto set X
%   of the objective functions defined in FITNESSFCN, subject to the linear
%   equalities Aeq*X = beq as well as the linear inequalities A*X <= b. (Set
%   A=[] and b=[] if no inequalities exist.) Linear constraints are supported
%   only for default PopulationType option ('doubleVector'); other population
%   types e.g., 'bitString' and 'custom' are not supported.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub) defines a set of
%   lower and upper bounds on the design variables, X, so that a local Pareto
%   set is found in the range lb <= X <= ub. Use empty matrices for lb and ub
%   if no bounds exist. Set lb(i) = -Inf if X(i) is unbounded below;  set
%   ub(i) = Inf if X(i) is unbounded above. Bound constraints are
%   supported only for default PopulationType option ('doubleVector');
%   other population types e.g., 'bitString' and 'custom' are not
%   supported.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,NONLCON) finds a local 
%   Pareto set X that satisfy the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. The Pareto set X must 
%   be such that C(X)<=0 and Ceq(X)=0. Set lb=[] and/or ub=[] if no bounds 
%   exist. Nonlinear constraints are not satisfied when the PopulationType 
%   option is set to 'bitString' or 'custom'. See the documentation for details.
%
%   X = GAMULTIOBJ(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,NONLCON,options)
%   finds a Pareto set X with the default optimization parameters replaced
%   by values in OPTIONS. OPTIONS can be created with the OPTIMOPTIONS
%   function. See OPTIMOPTIONS for details. For a list of options accepted
%   by GAMULTIOBJ refer to the documentation.
%
%   X = GAMULTIOBJ(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure that has the following fields:
%       fitnessfcn: <Fitness function>
%            nvars: <Number of design variables>
%            Aineq: <A matrix for inequality constraints>
%            bineq: <b vector for inequality constraints>
%              Aeq: <Aeq matrix for equality constraints>
%              beq: <beq vector for equality constraints>
%               lb: <Lower bound on X>
%               ub: <Upper bound on X>
%          nonlcon: <Nonlinear constraint function>
%          options: <Options created with optimoptions('gamultiobj',...)>
%           solver: <solver name 'gamultiobj'>
%         rngstate: <State of the random number generator>
%
%   [X,FVAL] = GAMULTIOBJ(FITNESSFCN,NVARS, ...) in addition returns a
%   matrix
%   FVAL, the value of all the objective functions defined in FITNESSFCN at
%   all the solutions in X. FVAL has numberOfObjectives columns and same number
%   of rows as does X.
%
%   [X,FVAL,EXITFLAG] = GAMULTIOBJ(FITNESSFCN,NVARS, ...) in addition returns
%   EXITFLAG which describes the exit condition of GAMULTIOBJ. Possible values
%   of EXITFLAG and the corresponding exit conditions are
%
%     1 Average change in value of the spread of Pareto set over
%        options.MaxStallGenerations generations less than options.FunctionTolerance.
%     0 Maximum number of generations exceeded.
%    -1 Optimization terminated by the output or plot function.
%    -2 No feasible point found.
%    -5 Time limit exceeded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = GAMULTIOBJ(FITNESSFCN,NVARS, ...) in addition
%   returns a structure OUTPUT with the following fields:
%            rngstate: <State of the random number generator before GA started>
%         generations: <Total number of generations, excluding HybridFcn iterations>
%           funccount: <Total number of function evaluations>
%       maxconstraint: <Maximum constraint violation>, if any
%             message: <GAMULTIOBJ termination message>
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION] = GAMULTIOBJ(FITNESSFCN, ...) in
%   addition returns the final POPULATION at termination. 
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION,SCORE] = GAMULTIOBJ(FITNESSFCN, ...)
%   in addition returns the SCORE of the final POPULATION.
%
%
%   Example:
%    Multiobjective minimization of two functions 'ackleyfcn' and 'shufcn'
%
%      fun1 = @(x) ackleyfcn(x);
%      fun2 = @(x) shufcn(x);
%      % Combine two objectives 'fun1' and 'fun2' 
%      fun1and2 = @(x) [fun1(x) fun2(x)];
%      % Bound constraints on X
%      lb = [-10 -10]; ub = [10 10];
%      % Specify the initial range for population
%      options = optimoptions('gamultiobj','InitialPopulationRange',[lb;ub]);
%      % Minimize using GAMULTIOBJ
%      [x,fval] = gamultiobj(fun1and2,2,[],[],[],[],lb,ub,options)
%
%    Display Pareto fronts and rank of individuals while GAMULTIOBJ
%    minimizes 
%
%      options = optimoptions('gamultiobj','InitialPopulationRange',[lb;ub]);
%      options = optimoptions(options,'PlotFcn',{@gaplotpareto,@gaplotrankhist});
%      [x,fval,exitflag,output] = gamultiobj(fun1and2,2,[],[],[],[],lb,ub,options)
%
%
%   See also OPTIMOPTIONS, PATTERNSEARCH, GA, @.

%   Reference: Kalyanmoy Deb, "Multi-Objective Optimization using
%   Evolutionary Algorithms", John Wiley & Sons ISBN 047187339


%   Copyright 2007-2015 The MathWorks, Inc.

% If the first arg is not a gaoptimset, then it's a fitness function followed by a genome
% length. Here we make a gaoptimset from the args.
defaultopt = struct('PopulationType', 'doubleVector', ...
    'PopInitRange', [], ...
    'PopulationSize', '50 when numberOfVariables <= 5, else 200', ...
    'CrossoverFraction', 0.8, ...
    'ParetoFraction', 0.35, ...
    'MigrationDirection','forward', ...
    'MigrationInterval',20, ...
    'MigrationFraction',0.2, ...
    'Generations', '200*numberOfVariables', ...
    'TimeLimit', inf, ...
    'StallGenLimit', 100, ...
    'TolFun', 1e-4, ...
    'TolCon', 1e-3, ...
    'InitialPopulation',[], ...
    'InitialScores', [], ...
    'PlotInterval',1, ...
    'CreationFcn',@gacreationuniform, ...
    'SelectionFcn', {{@selectiontournament,2}}, ...
    'CrossoverFcn',@crossoverintermediate, ...
    'MutationFcn',@mutationadaptfeasible, ...
    'DistanceMeasureFcn',{{@distancecrowding, 'phenotype'}}, ...
    'HybridFcn',[], ...
    'Display', 'final', ...
    'PlotFcns', [], ...
    'OutputFcns', [], ...
    'Vectorized', 'off', ...
    'UseParallel', false);

numOutputsRequested = nargout;

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && numOutputsRequested <= 1 && strcmpi(fun,'defaults')
    x = defaultopt;
    return
end

% Check number of input arguments
try 
    narginchk(1,10);
catch ME
    error(message('globaloptim:gamultiobj:numberOfInputs', ME.message));
end

if nargin < 10, options = [];
    if nargin < 9, nonlcon = [];
        if nargin < 8, ub = [];
            if nargin < 7, lb = [];
                if nargin <6, beq = [];
                    if nargin <5, Aeq = [];
                        if nargin < 4, bineq = [];
                            if nargin < 3, Aineq = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% For backwards compatibility, check to see if the options were passed
% in where "nonlcon" currently is.
if ~isempty(nonlcon) && (isstruct(nonlcon) || isa(nonlcon, 'optim.options.SolverOptions')) && nargin < 10
    options = nonlcon;
    nonlcon = [];
end

% One input argument is for problem structure
if nargin == 1
    if isa(fun,'struct')
        [fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,rngstate,options] = separateOptimStruct(fun);
        % Reset the random number generators
        resetDfltRng(rngstate);
    else % Single input and non-structure.
        error(message('globaloptim:gamultiobj:invalidStructInput'));
    end
end

% Prepare the options for the solver
options = prepareOptionsForSolver(options, 'gamultiobj');

% If fun is a cell array with additional arguments get the function handle
if iscell(fun)
    FitnessFcn = fun{1};
else
    FitnessFcn = fun;
end
% Only function handles or inlines are allowed for FitnessFcn
if isempty(FitnessFcn) ||  ~(isa(FitnessFcn,'inline') || isa(FitnessFcn,'function_handle'))
    error(message('globaloptim:gamultiobj:needFunctionHandle'));
end

% We need to check the nvars here before we call any solver
valid = isnumeric(nvars) && isscalar(nvars) && (nvars > 0) ...
    && (nvars == floor(nvars));
if ~valid
    error(message('globaloptim:gamultiobj:notValidNvars'));
end

user_options = options;
defaultopt.PopInitRange = [-10;10];
% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
        error(message('globaloptim:gamultiobj:optionsNotAStruct'));
elseif isempty(options)
    options = defaultopt;
end

% Take defaults for parameters that are not in options structure
options = gaoptimset(defaultopt,options);

% If a user doesn't specify PopInitRange, we want to set it to the
% bounds when we create the initial population. Need to store a flag
% that indicates whether the user has specified PopInitRange so we can
% do this in the creation function.
options.UserSpecPopInitRange = isa(user_options, 'struct') && ...
    isfield(user_options, 'PopInitRange') && ~isempty(user_options.PopInitRange);

% Check for non-double inputs
msg = isoptimargdbl('GAMULTIOBJ', {'NVARS','A',   'b',   'Aeq','beq','lb','ub'}, ...
                                    nvars,  Aineq, bineq, Aeq,  beq,  lb,  ub);
if ~isempty(msg)
    error('globaloptim:gamultiobj:dataType',msg);
end

% Introduce field to describe the objective type
options.MultiObjective = true;

% Validate constraints and options
[x,fval,exitFlag,output,population,scores,FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
    ~,options] = gacommon(nvars,fun,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,[],options,user_options);

if exitFlag < 0 % Infeasibility
    return;
end

if ~isempty(options.OutputFcns) || ~isempty(options.PlotFcns)
    % For calling an OutputFcn, make an options object (to be updated
    % later) that can be passed in
    options.OutputPlotFcnOptions = optimoptions(@gamultiobj);
    options.OutputPlotFcnOptions = copyForOutputAndPlotFcn(options.OutputPlotFcnOptions,options);
end

% Call appropriate single objective optimization solver
[x,history,fval,exitFlag,output,population,scores] = gamultiobjsolve2(FitnessFcn,nvars, ...
     Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options,output,numOutputsRequested);

 
