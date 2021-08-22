function Best = hsa1(CostF, Low, High, opts)
% Harmony Search Algorithm
%By Sajjad Yazdani
%
% Base on:
% [1]:
% Problem Prametters
fig = figure;

%stop button
pObj = uicontrol(fig,'string','Stop', ... 
      'Position',[2 5 50 20],'callback',{@Push});

han=gca;
% han=subplot(2,1,1);
% han2=subplot(2,1,2);

Dim=length(Low); % problem Dimention
Low=Low'; % Low Boundry of Problem
High=High'; % High Boundry of Problem
Min=1; % Minimaization or maximaiz of Fun? if Min=1 it will be minimaze the function and if Min=0 it will be maximized the function.
% Harmony Search Parametters
HMS=opts.population;
bw=opts.bw;
HMCR=opts.HMCR;
PAR=opts.PAR;
MaxItr=opts.MaxItr;

BestCost=zeros(MaxItr,1);    % Array to Hold Best Cost Values
Vars=BestCost;

BestSol.Cost=inf;
bestGlobal.Cost=inf;
bestGlobal.Time=inf;

% Initialization
HM=zeros(HMS,Dim);
HF=zeros(HMS,1);
for i=1:HMS
    HM(i,:)=Low+(High-Low).*rand(1,Dim);
    HF(i,1)=CostF(HM(i,:));
end
if Min==1
    [WorstFit,WorstLoc]=max(HF);
else
    [WorstFit,WorstLoc]=min(HF);
end

user_stopped=0;
% Iteration Loop
Itr=1;
while user_stopped==0 && Itr <= MaxItr
% for Itr=1:MaxItr
    
    i=1;
    while user_stopped==0 && i <= HMS
%     for i=1:HMS
    HarmonyIndex=fix(rand(1,Dim)*HMS)+1;% Random Selection of Harmony
    Harmony=diag(HM(HarmonyIndex, : ))';% Extraxt Value of harmony from Memory(Can Be better???)
    CMMask=rand(1,Dim)<HMCR;
    NHMask=(1-CMMask);
    PAMask=(rand(1,Dim)<PAR).*(CMMask);
    CMMask=CMMask.*(1-PAMask);
    NewHarmony=CMMask.*Harmony+PAMask.*( Harmony + bw*(2*rand(1,Dim)-1) )+NHMask.*(Low+(High-Low).*rand(1,Dim));
    OutOfBoundry=(NewHarmony>High)+(NewHarmony<Low);
    NewHarmony(OutOfBoundry==1)=Harmony(OutOfBoundry==1);
    NHF=CostF(NewHarmony);
    if (NHF<WorstFit)&&(Min==1)
        HM(WorstLoc,:)=NewHarmony;
        HF(WorstLoc)=NHF;
        [WorstFit,WorstLoc]=max(HF);
    elseif (NHF<WorstFit)&&(Min==0)
        HM(WorstLoc,:)=NewHarmony;
        HF(WorstLoc)=NHF;
        [WorstFit,WorstLoc]=min(HF);
    end
    i=i+1;
    end
    %{
    if Min==1
    [BestFit,BestLoc]=min(HF);
    else
    [BestFit,BestLoc]=max(HF);
    end
    %}
    [BestSol,BestLoc]=min(HF);
    
    if BestSol < bestGlobal.Cost
        bestGlobal.Cost = BestSol;
        bestGlobal.Time = toc;
        bestGlobal.Iter = Itr;
    end
    BestCost(Itr)=bestGlobal.Cost;
    %visualize costs
    drawnow

    plot(han,1:Itr,BestCost(1:Itr),'LineWidth',2);
    xlabel(han,'Iteration');
    ylabel(han,'Best Cost');
    title(han,sprintf('HSO iteration: %d',Itr))
    grid(han,'on');
%     han.YScale='log';

    
%     Vars(Itr) = sum(var(HM))/HMS;
%     plot(han2,1:Itr,Vars(1:Itr),'LineWidth',2);

    Itr=Itr+1;
end
% Present Best Answer
Best=HM(BestLoc,:);

function Push(source,event)
user_stopped=1;
end
end

