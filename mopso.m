function [REP,swarm]= mopso(objective,lower_bound,upper_bound,constraint,options)
%mopso is an implementation of multi objective particle swarm optimization
%technique for a minimization problem
%   when called mopso() it solves the provided example
%% initialize parameters
if nargin==0
    c = [0.1,0.2]; % [cognitive acceleration, social acceleration] coefficients
    iw = [0.5 0.001]; % [starting, ending] inertia weight
    max_iter = 200; % maximum iterations 100
    lower_bound = zeros(1,3); % lower bound of vars
    upper_bound = pi/2*ones(1,3); % upper bound of vars
    swarm_size=400; % swarm size 100
    rep_size=400; % Repository Size 100
    grid_size=70; % Number of Grids per Dimension 7
    alpha=0.1; % Inflation Rate
    beta=2; % Leader Selection Pressure
    gamma=2; % Deletion Selection Pressure
    mu=0.1; % Mutation Rate
    objective=@fitness; % objective function
    constraint=@constraints; % constraints function
else
    c = options.c;
    iw = options.iw;
    max_iter = options.max_iter;
    swarm_size = options.swarm_size;
    rep_size = options.rep_size;
    grid_size = options.grid_size;
    alpha = options.alpha;
    beta = options.beta;
    gamma = options.gamma;
    mu = options.mu;
end
%% initialize particles
fprintf('Initializing swarm ...\n')
w = @(it) ((max_iter - it) - (iw(1) - iw(2)))/max_iter + iw(2);
pm = @(it) (1-(it-1)/(max_iter-1))^(1/mu);
swarm(1,swarm_size) = Particle();
for i =1:swarm_size
    swarm(i)=Particle(lower_bound,upper_bound,objective,constraint);
    retry = 0;
    while ~all(swarm(i).isFeasable) && retry<100
        swarm(i)=Particle(lower_bound,upper_bound,objective,constraint);
        retry = retry + 1;
    end
end
REP = Repository(swarm,rep_size,grid_size,alpha,beta,gamma);
%% Loop
fprintf('Starting the optimization loop ...\n')
for it=1:max_iter
%         leader = REP.SelectLeader();
    %override leader finding, to let it use minimum
    costs=vertcat(swarm.cost); % get costs
    [~,leader]=min(costs);    
    leader=swarm(leader(1));
    
    wc = w(it); %current inertia weight
    pc=pm(it); %current mutation rate
    for i =1:swarm_size %update particles
        swarm(i)=swarm(i).update(wc,c,pc,leader,objective,constraint);
    end
%     feasable = all(horzcat(swarm.isFeasable)<0); % feasibility check!!
%     
%     % change infeasible to feasible
%     feas = sum(feasable);
%     
%     if feas < swarm_size
%         num = swarm_size-feas;
%         missing = find(~feasable);
%         good = find(feasable);
%         ids = randperm(feas);
%         for n = 1:num
%             % find missing adress and missing(n) to copy from good address
%         goodID = good(ids(mod(n-1,feas)+1)); % random using permutaion
%         swarm(missing(n)).x = swarm(goodID).x;
%         swarm(missing(n)).v = zeros(1,length(upper_bound));
%         feasable(missing(n)) = all(constraint(swarm(missing(n)).x)<0);
%         end
%     end
%     
%     if any(feasable)
%         last_feasible = swarm(feasable);
%     else
%         swarm = last_feasible;
%         REP = REP.update(swarm);
%         Title = sprintf('Iteration %d, Number of Rep Members = %d',it,length(swarm));
%         PlotCosts(swarm,swarm,Title)
%         fprintf('\nOptimisation stopped too early after %d iterations, no more feasible particles...\n',it)
%         return
%     end
     Title = sprintf('Iteration %d',it);
    PlotCosts(swarm,Title)
    %disp(Title);
end
end
function PlotCosts(swarm,Title)
figure(1)

feasable = horzcat(swarm.isFeasable)<0; % feasibility check!!
feasable_swarm = swarm(feasable);

LEG = {};
if ~isempty(feasable_swarm)
    swarm_costs=vertcat(feasable_swarm.cost);
    best_cost=min(swarm_costs);
    swarm_costs=swarm_costs(swarm_costs<best_cost*100); % take up to 100x min only
    n=numel(swarm_costs);
    plot(swarm_costs(1:n/2),swarm_costs(n/2+1:n),'go')
    hold on
    LEG = {'Current feasible SWARM'};
    Title = sprintf([Title '\nfeasable swarm=%d'],n/2);

%     Title = sprintf([Title '\nfeasable swarm=%d'],length(feasable_swarm));
end

xlabel('1^{st} Objective')
ylabel('2^{nd} Objective')
grid on
hold off
title(Title)
legend(LEG,'location','best')
drawnow
end
