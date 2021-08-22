%% Initialize variables
clear global
clear variables
close all
% warning('mopso.m does not have early stopping implemented yet')
% warning('BuildGroupFunction.m is not self compiling')


% load constants from excel:
% LOADSAVE
load('Excel_variables.mat')
load('UCs.mat')
load('UBs.mat')

% load('NumericalValues.mat');
tff=0;  L=0;    ggg=0;
Ned=0;  MYed=0; MZed=0; Ved=0;
moment1=0;  moment3=0;  moment5=0;  

%displacement values
indexes=1:GDof;
%can be generalized for N dimensions but whatever... we don't build tesseracts
r=[ [4:6:GDof] [5:6:GDof] [6:6:GDof]];
%removing rotations
indexes(r)=[];

%precompute DoF we are interested in
activeDof=(1:GDof)';
logUA = ~(ismember(activeDof,prescribedDof));
activeDof = activeDof(logUA);

% constants for visualization
Ncomb=9;	Nobjectives=2;
candidate=cell(Ncomb,1);    Plots.figs=cell(Ncomb,1);   Plots.names=cell(Ncomb,1);
times=zeros(Ncomb,1); bestTimes=zeros(Ncomb,2,Nobjectives);   costHistory=cell(Ncomb,2);    figNum=0;
                     %^ 1st-time, 2nd-iter ^

% data to send into functions
params = struct('activeDof',activeDof,'beta',beta ,'C1',C1 ,'DisloadDirection',DisloadDirection,'E',E,...
    'elementDof',elementDof,'epsilon',epsilon,'forceVector',forceVector,'G',G,...
    'elementNodes',elementNodes,'ggg',ggg,'gammaM0',gammaM0,'gammaM1',gammaM1,...
    'indexes',indexes, 'ita',ita,  'KgElement',0,'Ved',Ved, 'MYed',MYed,...
    'MZed',MZed, 'll',ll,'LiveL',0, 'Load1',Load1, 'loadpoints',loadpoints,'L',L,...
    'lambdaDashLT0',lambdaDashLT0   ,'moment1',moment1 ,'moment3',moment3,...
    'moment5',moment5,'Ned',Ned,'phiColumnYY',0  ,'phiColumnZZ',0  ,'phiLT',0,...
    'psai',0 ,'psi',0  ,'ki',0,'memberlength',memberlength,'nodeCoordinates',nodeCoordinates,...
    'numberElements',numberElements, 'stiffness',0, 'T',0, 'Tt',0,'loadmembers',loadmembers,...
    'loadjoints',0,'GDof',GDof,'steelGrade',steelGrade,...
    'V',V,'xx',xx,'yy',yy,'zz',zz,'L0',L0,'UCs',UCs,'UBs',UBs);

afterParams = struct('continuity',0, 'grouping',1, 'is2D',0, 'multiObj',0, 'Nobjectives',Nobjectives, 'type',[]);

%% Optimisation part


A=[1 5 9 13 4 8 12 16 17 21 25 29 20 24 28 32 33 37 41 45 36 40 44 48];
B=[121 125 129 133 124 128 132 136 137 141 145 149 140 144 148 152 153 157 161 165 156 160 164 168];    %6STOREY
C=[2 6 10 14 3 7 11 15 18 22 26 30 19 23 27 31 34 38 42 46 35 39 43 47];
D=[122 126 130 134 123 127 131 135 138 142 146 150 139 143 147 151 154 158 162 166 155 159 163 167];
E=[49:120 169:216];
F=[217:240];

Order={A,B,C,D,E,F};


BuildGroupFunction(Order)
nvars=length(Order); % number of groups
intcon = 1:nvars;

ub=zeros(nvars,1);
for id=1:nvars %read bounds from corresponding groups
   ub(id) =  46*(SectionType(Order{id}(1)) == 2) + 107 * (SectionType(Order{id}(1)) == 1);
end
lb = ones(nvars,1);

% Single Objective definitions
singleoptions=gaoptimset('Generations',100,'EliteCount',100,'SelectionFcn',@selectiontournament, ...
    'PopulationSize',200,'CrossoverFraction',0.7,'FitnessScalingFcn', @fitscalingrank, ...
    'PlotFcns', {@gaplotbestf, @gaplotbestindiv, @gaplotexpectation});

psooptions = optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'Display','final','Plotfcn',@pswplotbestf);

hsoptions = struct('population', 200, 'bw', 5, 'HMCR', 0.6,'PAR', 0.5, 'MaxItr', 100);


% Multi Objective definitions
% {@gaplotbestf, @gaplotbestindiv, @gaplotexpectation,@gaplotdistance,@gaplotgenealogy, @gaplotselection , @gaplotscorediversity, @gaplotrange});

% !!! renamed gaplotpareto function so that it works as we want !!!
multioptions = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto2,'Generations',100,'TolCon',1e-7,'SelectionFcn',@selectiontournament, ...
    'CrossoverFcn',@crossoverheuristic,'PopulationSize',100,'HybridFcn',@fgoalattain,'Display','final');

mopsoptions = struct('c', [0.1,0.2], 'iw', [0.5 0.001], 'max_iter', 100, 'swarm_size', 100, 'grid_size', 70,...
    'alpha', 0.1, 'beta', 2 , 'gamma', 2, 'mu', 0.1);
mohsoptions = struct('population', 100, 'bw', 0.3, 'HMCR', 0.3,'PAR', 0.3, 'MaxItr', 100);

% User interaction
clc
Choice = input(['Which parameters optimize to?\n 1 - weight\n 2 - energy\n 3 - multi (weight+energy)\n'...
    'H - HSA, E - mod HSA,  G - Ga, P - PSO\nExample: 1G - GA weight, 3P - PSO mixed,  2EGHP-runs every model on energy cost...\n>> '],'s');

algoT = {'PSO','GA','HSA','EGHS'};
goalT = {'weight','energy','Multi-obj'};
goalO = 'single';

algoID=find( ismember('PGHE',Choice) );
%goalID=find( ismember('123',Choice) );
goalID=Choice(1)-48;

nAlgs=length(algoID);
nGoal=length(goalID);

for iGoal=goalID % for future
    for iAlg=algoID
        
        if iGoal==3,    goalO = 'multi';    end
        
        fprintf('\nStarting %s objective %s - %s...\n',goalO,algoT{iAlg},goalT{iGoal})
        tic
        
        ConF=@(x) constraints(group(round(x)),params);
        
        if iGoal==3
            
            CostMOPSO=@(row)myCost(UCs,UBs,group(round(row)),params,103);
            CostGAMulti=@(row)myCost(UCs,UBs,group(round(row)),params,103);
            CostMOHSA=CostMOPSO;

			afterParams.multiObj=1;
            
            switch iAlg
                case 4
                    error('MOHSA integer not implemented yet')
                case 3
                    [x_best,history] = mohsa(CostMOHSA,lb,ub,mohsoptions);
                case 2
                    [x_best,history] = gamultiobj2(CostGAMulti,nvars,[],[],[],[],lb,ub,ConF,multioptions);
                                        
                    % Find unique from population
                    x_gaMul = unique(round(x_best),'rows');
                    
                    x_best=x_gaMul(1,:); %group done in afterOptim
                case 1
%                     ConF=@(x) -1; %workaround in cost function, since the unfeasible sets are ignored
                    % bounds transpositions due to way mopso was coded
                    
%                     [~,Swarm,history] = mopso(CostMOPSO,lb',ub',ConF,mopsoptions);
                    [best,history,swarm] = mopso(CostMOPSO,lb',ub',ConF,mopsoptions);

                     x_best = round(best.x);
%                     x_best = round(Swarm(1).x);
            end
        else
            CostGA=@(row)myCost(UCs,UBs,group(row),params,iGoal); %last param either 1 or 2
            CostPSO=@(row)myCost(UCs,UBs,group(round(row)),params,100+iGoal);
            CostHSA=CostPSO;
			afterParams.multiObj=0;
            
            switch iAlg %call algorithm
                case 4
                    x_best = hsa_integer(CostHSA,lb,ub,hsoptions);
                    
                case 3
                    x_best = hsa1(CostHSA,lb,ub,hsoptions);
                    
                case 2
                    x_best = ga(CostGA,nvars,[],[],[],[],lb,ub,ConF,intcon,singleoptions);
                    
                    le=length(singleoptions.PlotFcns); % check which subplot we show
                    subplot(le,1,1)
                    
                case 1
%                                         x_best =
%                                         aco2(CostACO,lb,ub,acoptions); %
%                                         for comparison with deleted index

                    x_best = round(particleswarm(CostPSO,nvars,lb,ub,psooptions));
            end
        end%if iGoal==3
		
        afterParams.type=sprintf( '%s %s', algoT{iAlg}, goalT{iGoal} );
        afterOptim(afterParams)
    end%for iAlg
end%for iGoal
%% Visualization part
Range = 1:figNum;
visualize(figNum,Plots,cellfun(group,candidate(Range),'UniformOutput',0),params,afterParams.continuity)

%% Comparison
comparison(goalID,costHistory,bestTimes,Plots,Range)
%% Helpful functions
%moved to separate files