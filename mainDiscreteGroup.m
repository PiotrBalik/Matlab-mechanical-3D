%% Initialize variables
clear global
clear variables
close all

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
Ncomb=9;
candidate=cell(Ncomb,1);    Plots.figs=cell(Ncomb,1);   Plots.names=cell(Ncomb,1);
times=cell(Ncomb,1);    figNum=0;


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

%% Optimisation part

A=[1 2 3 4];
B=[5 6 7 8];

% A=[1:8 19:24 32:33];
% B=setdiff(1:36,A);
% C=[15 21 22 28];
% D=[17 19 24 26];
% E=[29 35 36 42];
% F=[31 33 38 40];
% G=[43 49 50 56];
% H=[45 47 52 54];
% I=[2 4 6 9 11 13 16 18 20 23 25 27 30 32 34 37 39 41 44 46 48];
% J=[51 53 55];

%ordered indexes for easier codes :)
% Order={A, B, C, D, E, F, G, H, I, J};
Order={A,B};

BuildGroupFunction(Order)
nvars=2; % manual number of group
intcon = 1:nvars;
% ub = 46*(SectionType == 2) + 107 * (SectionType == 1);
ub=[46 107]'; % each group upper bounds
lb = ones(nvars,1);

% Single Objective definitions

singleoptions=gaoptimset('Generations',300,'EliteCount',100,'SelectionFcn',@selectiontournament, ...
    'PopulationSize',200,'CrossoverFraction',0.7,'FitnessScalingFcn', @fitscalingrank, ...
    'PlotFcns', {@gaplotbestf, @gaplotbestindiv, @gaplotexpectation});

psooptions = optimoptions('particleswarm','SwarmSize',250,'HybridFcn',@fmincon,'Display','final','Plotfcn',@pswplotbestf);

acoptions = struct('MaxIt',100, 'nAnt',150, 'Q',2,'tau0', 100/nvars , 'alpha',2, 'beta', 1, 'rho', 0.02);


ConGA=@(x) constraints(group(x),params); %define constraints function for GA

CostGA_M=@(row)myCost(UCs,UBs,group(row),params,1);

CostGA_E=@(row)myCost(UCs,UBs,group(row),params,2);

CostPSO_M=@(row)myCost(UCs,UBs,round(group(row)),params,101); %ACO uses PSO, since they handle constraints the same

CostPSO_E=@(row)myCost(UCs,UBs,round(group(row)),params,102);


% Multi Objective definitions

multioptions = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto,'Generations',300,'TolCon',1e-7,'SelectionFcn',@selectiontournament, ...
    'CrossoverFcn',@crossoverheuristic,'PopulationSize',200,'HybridFcn',@fgoalattain,'Display','final');

mopsoptions = struct('c', [0.1,0.2], 'iw', [0.5 0.001], 'max_iter', 200, 'swarm_size', 400, 'grid_size', 70,...
                        'alpha', 0.1, 'beta', 2 , 'gamma', 2, 'mu', 0.1);

ConGAMulti=@(x) constraints(round(group(x)),params); %define constraints function for GA multi

CostGAMulti=@(row)myCost(UCs,UBs,round(group(row)),params,103); %define cost function fo GA multi

CostMOPSO=@(row)myCost(UCs,UBs,round(group(row)),params,103);

% {@gaplotbestf, @gaplotbestindiv, @gaplotexpectation,@gaplotdistance,@gaplotgenealogy, @gaplotselection , @gaplotscorediversity, @gaplotrange});

clc
 Choice = input(['Which parameters optimize to?\n 1 - weight\n 2 - energy\n 3 - multi (weight+energy)\n'...
     'A - ACO, G - Ga, P - PSO\nExample: 1G, 31P,  123GAP-runs every model...\n>> '],'s');

if any(Choice=='1')
    if any(Choice=='A')
        
        disp('Starting single objective ACO - mass...')
        tic
        
        x_best = aco(CostPSO_M,lb,ub,acoptions);
        
        afterOptim('ACO weight')
        
        x_acoM=x_best;
    end
    if any(Choice=='G')
        
        disp('Starting single objective GA  - mass...')
        tic; % start real-time measure
        
        x_best=ga(CostGA_M,nvars,[],[],[],[],lb,ub,ConGA,intcon,singleoptions);
        
        le=length(singleoptions.PlotFcns); % check which subplot we show, to title it
        subplot(le,1,1)
        afterOptim('GA weight')
    end
    if any(Choice=='P')
        
        disp('Starting single objective PSO - mass...')
        tic
        
        x_best = round(particleswarm(CostPSO_M,nvars,lb,ub,psooptions));
        
        afterOptim('PSO weight')
        x_psoM=x_best;
    end
end%choice==1

if any(Choice=='2')
     if any(Choice=='A')
        
        disp('Starting single objective ACO - energy...')
        tic
        
        x_best = aco(CostPSO_E,lb,ub,acoptions);
        
        afterOptim('ACO energy')
        
        x_acoE=x_best;
    end
    if any(Choice=='G')
        
        disp('Starting single objective GA - energy...')
        tic
        
        x_best=ga(CostGA_E,nvars,[],[],[],[],lb,ub,ConGA,intcon,singleoptions);
        
        le=length(singleoptions.PlotFcns); % check which subplot we show
        subplot(le,1,1)
        afterOptim('GA energy')
        x_gaE = x_best;
    end
    if any(Choice=='P')
        
        disp('Starting single objective PSO - energy...')
        tic;
        
        x_best = round(particleswarm(CostPSO_E,nvars,lb,ub,psooptions));
        
        afterOptim('PSO energy')
        x_psoE=x_best;
    end
    
end% choice==2

if any(Choice=='3')
    if any(Choice=='G')
        
        disp('Starting Multi objective GA...')
        tic
        
        x_best=gamultiobj(CostGAMulti,nvars,[],[],[],[],lb,ub,ConGAMulti,multioptions);
        
        x_gaMul=round(x_best);
        
        % Remove unique from founding
        x_gaMul = unique(x_gaMul,'rows');
        
        %{
        [m,~]=size(x_gaMul);
        for a=1:m % print every unique candidate
            myCost(UCs,UBs,group(x_gaMul(a,:)),params,10);
        end
        %}
        x_best=x_gaMul(1,:); %group done in afterOptim
        afterOptim('GA multi')
    end
    if any(Choice=='P')
        
        disp('Starting Multi objective PSO...')
        tic
        
        ConF=@(x) -1; %workaround in cost function, since the unfeasible sets are ignored
        % bounds transpositions due to way mopso was coded
        [Rep,Swarm] = mopso(CostMOPSO,lb',ub',ConF,mopsoptions);
        
        x_best = round(Swarm(1).x);
        afterOptim('PSO multi')
    end
    
end% choice==3

%% Visualization part
continuity=0;
visualize(figNum,Plots,candidate,params,continuity)

%% Comparison
%candidates: for id=1:figNum, disp(candidate{id}); disp(times{id}) end
%table of times:

%% Helpful functions
%moved to files