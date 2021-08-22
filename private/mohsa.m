function [Best, history] = mohsa(CostF, Low, High, opts)
% Harmony Search Algorithm
%By Sajjad Yazdani
%
% Base on:
% [1]:
% Problem Prametters
figure;
han=gca;

Dim=length(Low); % problem Dimention
Low=Low'; % Low Boundry of Problem
High=High'; % High Boundry of Problem
% Harmony Search Parametters
HMS=opts.population;
bw=opts.bw;
HMCR=opts.HMCR;
PAR=opts.PAR;
MaxItr=opts.MaxItr;

history.Costs = zeros(MaxItr,2);
history.bestIt = zeros(1,2);
history.bestC = inf(1,2);
history.lastIt = 0;
history.bestTime = zeros(1,2);

% Initialization
HM=zeros(HMS,Dim);
HF=zeros(HMS,2);
for i=1:HMS
    HM(i,:)=Low+(High-Low).*rand(1,Dim);
    HF(i,:)=CostF(HM(i,:));
end
%worst pareto elements
[WorstFit,WorstLoc]=find_pareto_frontier(1./HF);
% Iteration Loop
for Itr=1:MaxItr
    
    for i=1:HMS
        HarmonyIndex=fix(rand(1,Dim)*HMS)+1;% Random Selection of Harmony
        Harmony=diag(HM(HarmonyIndex,1:Dim))';% Extraxt Value of harmony from Memory(Can Be better???)
        CMMask=rand(1,Dim)<HMCR;
        NHMask=(1-CMMask);
        PAMask=(rand(1,Dim)<PAR).*(CMMask);
        CMMask=CMMask.*(1-PAMask);
        NewHarmony=CMMask.*Harmony+PAMask.*( Harmony + bw*(2*rand(1,Dim)-1) )+NHMask.*(Low+(High-Low).*rand(1,Dim));
        OutOfBoundry=(NewHarmony>High)+(NewHarmony<Low);
        %clamp the outliers
        NewHarmony(OutOfBoundry==1) = Harmony(OutOfBoundry==1);
        NHF=CostF(NewHarmony);
        
        improved = NHF < 1./WorstFit; %from inverse pareto
        if any(improved)
            impBoth = improved(:,1) & improved(:,2);
            
            if any(impBoth)
                bestID = find(impBoth);
                Nboth = numel(bestID);
                bestID = bestID(randperm( Nboth, 1)); %get one from the improvement
            else
                improved = improved(:,1) | improved(:,2);
                bestID = find(improved);   %which worst are dominated
                Nimpr = numel(bestID);
                bestID = bestID(randperm( Nimpr, 1));
            end
            
            toSwap = find(WorstLoc);
            
            HM( toSwap(bestID), :) = NewHarmony;
            HF( toSwap(bestID), :) = NHF;
            [WorstFit,WorstLoc]=find_pareto_frontier(1./HF);
        end

    end
    
    
    for nOb=1:2 %update history
        if Itr==1
           history.Costs(1,nOb)=min(HF(:,nOb));
           history.bestC(nOb)=history.Costs(1,nOb);
        end
        
        if min(HF(:,nOb)) < history.bestC(nOb)
            history.bestC(nOb) = min(HF(:,nOb));
            history.bestIt(nOb)=Itr;
            history.bestTime(nOb)=toc;
        end
        history.Costs(Itr,nOb)=history.bestC(nOb);
    end
    history.lastIt=Itr;

    %visualize costs
    drawnow
    
    plot(han,HF(:,1),HF(:,2),'pm');
    
    xlabel(han,'Weight');
    ylabel(han,'Energy');
    title(han,sprintf('HSO iteration: %d',Itr))
    grid(han,'on');
    
end
[~,BestLoc]=find_pareto_frontier(HF);

% Present Best Answer
Best=HM(find(BestLoc,1),:);
end