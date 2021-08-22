function [best,history,swarm]= mopso(objective,lower_bound,upper_bound,constraint,options)
%mopso is an implementation of multi objective particle swarm optimization
%technique for a minimization problem
%   when called mopso() it solves the provided example
% initialize parameters

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
history.Costs = zeros(max_iter,2);
history.bestIt = zeros(1,2);
history.bestC = inf(1,2);
history.lastIt = 0;
history.bestTime = zeros(1,2);

best=Particle;
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
        swarm(i).isFeasable=all(constraint(swarm(i).x )< -1e-6);
    end
    REP.update(swarm);
%     [best,history]=PlotCosts(swarm,it,history);
     [best,history]=PlotCosts(swarm,it,history,best);

end

function [best,history]=PlotCosts(swarm,it,history,best)

drawnow
feasable = horzcat(swarm.isFeasable)==1;
feasable_swarm = swarm(feasable);

Title = sprintf('Iteration %d',it);

LEG = {};
if ~isempty(feasable_swarm)
    swarm_costs=vertcat(feasable_swarm.cost);
    
    for nOb=1:2 %update history
        if it==1
           history.Costs(it,nOb)=min(swarm_costs(:,nOb));
           history.bestC(nOb)=history.Costs(it,nOb);
        end
        
        if min(swarm_costs(:,nOb)) < history.bestC(nOb)
            [ history.bestC(nOb),id ] = min(swarm_costs(:,nOb));
            history.bestIt(nOb)=it;
            history.bestTime(nOb)=toc;
            best=feasable_swarm(id);
        end
        history.Costs(it,nOb)=history.bestC(nOb);
    end
    history.lastIt=it;
    
    n=numel(swarm_costs);
    
    %determine along 1st objective
    %     best_cost=sort(swarm_costs(:,1));
    %     best_cost=mean(best_cost(1:round(n/10)));
    
    %     swarm_costs=swarm_costs(swarm_costs(:,1)<best_cost*100,:); % take up to 100x min only
    plot(swarm_costs(:,1),swarm_costs(:,2),'go')
    hold on
    
    LEG = {'Current feasible SWARM'};
end
Title = sprintf([Title '\nfeasible swarm=%d'],n/2);
xlabel('1^{st} Objective')
ylabel('2^{nd} Objective')
grid on
hold off
title(Title)
legend(LEG,'location','best')

end
end