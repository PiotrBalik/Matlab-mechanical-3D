%% Initialize variables
clear
close all

%to call from command window after optimization
%>>seeForces(x_best,params)

% load constants from excel:
% LOADSAVE
load('Excel_variables.mat')
load('UCs.mat')
load('UBs.mat')

% load('NumericalValues.mat');
tff=0;
Ned=0;
MYed=0;
MZed=0;
Ved=0;
moment1=0;
moment3=0;
moment5=0;
L=0;
ggg=0;

%displacement values
indexes=1:144;

%removing rotations
r=mod(indexes,[4 5 6]')==0;
%can be generalized for N dimensions but whatever...
r1=indexes(r(1,:));
r2=indexes(r(2,:));
r3=indexes(r(3,:));
r=unique([r1 r2 r3]);
indexes(r)=[];

% constants for visualize
candidate=cell(6,1);
Plots.figs=cell(6,1);
Plots.names=cell(6,1);
figNum=0;


% data to send into functions
params = struct('beta',beta ,'C1',C1 ,'DisloadDirection',DisloadDirection,'E',E,...
    'elementDof',elementDof,'epsilon',epsilon,'forceVector',forceVector,'G',G,...
    'elementNodes',elementNodes,'ggg',ggg,'gammaM0',gammaM0,'gammaM1',gammaM1,...
    'indexes',indexes, 'ita',ita,  'KgElement',0,'Ved',Ved, 'MYed',MYed,...
    'MZed',MZed, 'll',ll,'LiveL',0, 'Load1',Load1, 'loadpoints',loadpoints,'L',L,...
    'lambdaDashLT0',lambdaDashLT0   ,'moment1',moment1 ,'moment3',moment3,...
    'moment5',moment5,'Ned',Ned,'phiColumnYY',0  ,'phiColumnZZ',0  ,'phiLT',0,...
    'psai',0 ,'psi',0  ,'ki',0,'memberlength',memberlength,'nodeCoordinates',nodeCoordinates,...
    'numberElements',numberElements, 'stiffness',0, 'T',0, 'Tt',0,'loadmembers',loadmembers,...
    'loadjoints',0,'GDof',GDof,'prescribedDof',prescribedDof,'steelGrade',steelGrade,...
    'V',V,'xx',xx,'yy',yy,'zz',zz,'L0',L0,'UCs',UCs,'UBs',UBs);

%% Optimisation part

A=[1:8 19:24 32:33];
B=setdiff(1:36,A);
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

ConGA=@(x) constraints(group(x),params); %define constraints function for GA

CostGA_M=@(row)myCost(UCs,UBs,group(row),params,1);

CostGA_E=@(row)myCost(UCs,UBs,group(row),params,2);

CostPSO_M=@(row)myCostPSOdiscrete(UCs,UBs,round(group(row)),params,1);

CostPSO_E=@(row)myCostPSOdiscrete(UCs,UBs,round(group(row)),params,2);


% Multi Objective definitions

multioptions = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto,'Generations',300,'TolCon',1e-7,'SelectionFcn',@selectiontournament, ...
    'CrossoverFcn',@crossoverheuristic,'PopulationSize',200,'HybridFcn',@fgoalattain,'Display','final');

ConGAMulti=@(x) constraints(round(group(x)),params); %define constraints function for GA multi

CostGAMulti=@(row)myCost(UCs,UBs,round(group(row)),params,3); %define cost function fo GA multi

CostMOPSO=@(row)myCostPSOdiscrete(UCs,UBs,round(group(row)),params,3);

% {@gaplotbestf, @gaplotbestindiv, @gaplotexpectation,@gaplotdistance,@gaplotgenealogy, @gaplotselection , @gaplotscorediversity, @gaplotrange});
% each element is in range 1-107 or 1-46
% so for Uc we have 1<= Uc <= 46
% for UB we have 1 <= UB <=107
clc
Choice = input('Which parameters optimize to?\n 1 - weight\n 2 - energy\n 3 - multi (weight+energy)\nG - Ga, P - PSO\nExample: 1G, 31P,  123GP-runs every model...\n','s');
% nvars = numberElements;
% nvars=10; % manual number of group
% intcon = 1:nvars;
% % ub = 46*(SectionType == 2) + 107 * (SectionType == 1);
% ub=[46 46 46 46 46 46 46 46 107 107]'; % each group upper bounds
% lb = ones(nvars,1);


if any(Choice=='1')
 if any(Choice=='G')
    
tic; % start real-time measure
disp('Starting single objective GA  - mass...')

x_best=ga(CostGA_M,nvars,[],[],[],[],lb,ub,ConGA,intcon,singleoptions);

le=length(singleoptions.PlotFcns); % check which subplot we show, to title it
subplot(le,1,1)
afterOptim('GA weight')

 end
 if any(Choice=='P')

disp('Starting single objective PSO - mass...')

tic;

x_best = particleswarm(CostPSO_M,nvars,lb,ub,psooptions);

afterOptim('PSO weight')

x_psoM=x_best;
 end
end

if any(Choice=='2')
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


x_best = particleswarm(CostPSO_E,nvars,lb,ub,psooptions);

afterOptim('PSO energy')
x_psoE=x_best;

     end

end

if any(Choice=='3')
     if any(Choice=='G')
    
disp('Starting Multi objective GA...')

tic

x_best=gamultiobj(CostGAMulti,nvars,[],[],[],[],lb,ub,ConGAMulti,multioptions);

x_gaMul=round(x_best);

% Remove unique from founding
x_gaMul = unique(x_gaMul,'rows');

[m,~]=size(x_gaMul);
for a=1:m % print every unique candidate
    myCost(UCs,UBs,group(x_gaMul(a,:)),params,10);
end
    x_best=x_gaMul(1,:); %group done in afterOptim
    afterOptim('GA multi')
    
     end
     if any(Choice=='P')
	 
options.c=[0.1,0.2]; % [cognitive acceleration, social acceleration] coefficients 0.1 0.2
options.iw = [0.5 0.001]; % [starting, ending] inertia weight 0.5 0.001
options.max_iter = 200; % maximum iterations 200
options.swarm_size = 400; % swarm size 400
options.rep_size = 400; % Repository Size 400
options.grid_size = 70;% Number of Grids per Dimension 70
options.alpha = 0.1; % Inflation Rate 0.1 
options.beta = 2; % Leader Selection Pressure 2
options.gamma = 2; % Deletion Selection Pressure 2
options.mu = 0.1; % Mutation Rate 0.1


tic

ConF=@(x) -1;
% bounds transpositions due to way mopso was coded
[Rep,Swarm] = mopso(CostMOPSO,lb',ub',ConF,options);

% Swarm = mopso(@(row)myCost(UCs,UBs,row,params,3),lb',ub',@(x) constraints(x,params),options);
% % Swarm = mopso(@(row)myCostContinous(UCCS,UBBS,row,params,3),lb',ub',@(x) constraintsContinuous(x,params),options)
% x_best = Swarm.Swarm(1).x;
x_best =Swarm(1).x;
afterOptim('PSO multi')
     end

end

%% Visualization part
continuity=0;
visualize(figNum,Plots,candidate,params,continuity)

%% Helpful functions

function printTime
t_m = toc; % measure time
fprintf('\nTime elapsed: %2.0fm %2.0fs\n',floor(t_m/60),t_m-60*floor(t_m/60));
end

function  afterOptim(type)
grid on
grid minor;
evalin('base','figNum=figNum+1;');
%save the figure for future
evalin('base','Plots.figs{figNum} = print(''-RGBImage'');');
evalin('base',['Plots.names{figNum} =''' type ''';']);
disp(['Done ' type])
printTime

%unwrap group to compute cost of candidate in main workspace, rounding does not make a difference for intGA, but called anyway to clean up codelines
evalin('base','x_best=group(round(x_best));')

%show candidate
% evalin('base','disp(x_best)');
evalin('base','myCost(UCs,UBs,x_best,params,10);')
evalin('base','candidate{figNum}=x_best;')
evalin('base','constr_values=constraints(x_best,params)+1; disp(''Constraint values:''); disp(reshape(constr_values,8,numberElements))')
end