function [Abest, bestHist, bestGlobal]=aco(CostF,lb,ub,options)
% ACO main code
% Interfaced for chosing an index of the table. Cost function is determined by the elements at i-th position from table.
% heavily modified - BASED ON:
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YOEA103
% Project Title: Ant Colony Optimization for Traveling Salesman Problem
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com

% Results
figure;

% Parameters
MaxIt=options.MaxIt;    % Maximum Number of Iterations
nAnt=options.nAnt;      % Number of Ants (Population Size)
Q=options.Q;            % Phormone placement strength
tau0=options.tau0;      % Initial Phromone
alpha=options.alpha;    % Phromone Exponential Weight
beta=options.beta;      % Heuristic Exponential Weight
rho=options.rho;        % Evaporation Rate

nVar=max(length(lb)-1,1);

%function wrapper
%j=RouletteWheelSelection(P)
Roulette=@(P) find(rand<=cumsum(P),1,'first');

% Initialization
nNodes=max(ub)-min(lb) + 1;
tau=tau0*ones(  nNodes,nNodes, nVar);   % Phormone Matrix, 3D since each pair of groups adds another

eta=1./tau;
warning('Function ACO.m has not yet "eta" parameter implemented');

BestCost=zeros(MaxIt,1);    % Array to Hold Best Cost Values

% Empty Ant
empty_ant.Tour=[];  empty_ant.Cost=[];

% Ant Colony Matrix
ant=repmat(empty_ant,nAnt,1);

% Best Ant
BestSol.Cost=inf;
bestGlobal.Cost=inf;
Abest=[];

% ACO Main Loop
for it=1:MaxIt
    
    % Move Ants
    for k=1:nAnt
        
        ant(k).Tour=randi([lb(1) ub(1)]);
        
        for id=1:nVar
            
            i=ant(k).Tour(end);
            
            P=tau(i,:,id).^alpha.*eta(i,:,id).^beta;
            
            P(ant(k).Tour)=0;
            P=P/sum(P);
            
            %clamping should never be activated
            j=min(ub(id+1), max(lb(id+1), Roulette(P) ));
            
            ant(k).Tour=[ant(k).Tour j];
        end
        
        ant(k).Cost=CostF(ant(k).Tour);
    end%for k
    
    % Update Phromones
    for k=1:nAnt
        
        tour=ant(k).Tour;
        
        tour=[tour tour(1)]; %#ok
        
        for id=1:nVar
            
            i=tour(id);
            j=tour(id+1);
            
            tau(i,j,id)=tau(i,j,id)+Q/ant(k).Cost;
        end
    end% for k
    
    % Evaporation
    tau=(1-rho)*tau;

	
	% Store Best Cost
    [BestSol,NAbest]=min([ant(:).Cost]);
    if BestSol < bestGlobal.Cost
        Abest = ant(NAbest).Tour;
        bestGlobal.Cost = BestSol;
        bestGlobal.Iter = it;
    end
    BestCost(it)=BestSol;
    
    drawnow
    subplot(2,1,1)
    pcolor(tau(:,:,1))
    title('\tau values in the domain - 1st layer')
    
    subplot(2,1,2)
    plot(BestCost(1:it),'LineWidth',2);
    xlabel('Iteration');
    ylabel('Best Cost');
    title(sprintf('ACO iteration: %d',it))
    grid on;
end% for maxit
bestHist = BestCost;
