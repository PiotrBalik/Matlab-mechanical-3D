function [REP,swarm]= mopso(objective,lower_bound,upper_bound,constraint,options)
%mopso is an implementation of multi objective particle swarm optimization
%technique for a minimization problem
%   when called mopso() it solves the provided example
% initialize parameters
warning('mopso.m has not yet implemented early stopping')
if nargin==0
    warning('Called without arguments')
    return
else
    c = options.c;
    iw = options.iw;
    max_iter = options.max_iter;
    swarm_size = options.swarm_size;
    grid_size = options.grid_size;
    alpha = options.alpha;
    beta = options.beta;
    gamma = options.gamma;
    mu = options.mu;
end
% initialize particles
fprintf('Initializing swarm ...\n')
figure
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
REP = Repository(swarm,grid_size,alpha,beta,gamma);
% Loop
fprintf('Starting the optimization loop ...\n')
for it=1:max_iter
    leader = REP.SelectLeader(); %hopefully fixed
    %override leader finding, to let it use minimum
%     costs=vertcat(swarm.cost); % get costs
%     [~,leader]=min(costs);    
%     leader=swarm(leader(1));
    
    wc = w(it); %current inertia weight
    pc=pm(it); %current mutation rate
    for i =1:swarm_size %update particles
        swarm(i)=swarm(i).update(wc,c,pc,leader,objective,constraint);
    end
    REP.update(swarm);
    Title = sprintf('Iteration %d',it);
    PlotCosts(swarm,Title)
end
end
function PlotCosts(swarm,Title)

drawnow
feasable = horzcat(swarm.isFeasable)<0; % feasibility check!!
feasable_swarm = swarm(feasable);

LEG = {};
if ~isempty(feasable_swarm)
    swarm_costs=vertcat(feasable_swarm.cost);
    n=numel(swarm_costs);
    
    %determine along 1st objective
    best_cost=sort(swarm_costs(:,1));
    best_cost=mean(best_cost(1:round(n/10)));
    
    swarm_costs=swarm_costs(swarm_costs(:,1)<best_cost*100,:); % take up to 100x min only
    plot(swarm_costs(:,1),swarm_costs(:,2),'go')
    hold on
    LEG = {'Current feasible SWARM'};
    Title = sprintf([Title '\nfeasable swarm=%d'],n/2);
end

xlabel('1^{st} Objective')
ylabel('2^{nd} Objective')
grid on
hold off
title(Title)
legend(LEG,'location','best')

end
