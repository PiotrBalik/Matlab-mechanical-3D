function comparison(goalID,costHistory,bestTimes,Plots,Range)
figure
if goalID==3 %multiobjective comparison
    maxX = cellfun(@max,costHistory(Range,1));
    maxY = cellfun(@max,costHistory(Range,2));
    minX = cellfun(@min,costHistory(Range,1));
    minY = cellfun(@min,costHistory(Range,2));
    dX=max(maxX-minX);
    dY=max(maxY-minY);
    
    Cols=hsv(Range(end));
        for id=Range
        Leg{id}=plot(rand(2,1),-rand(2,1),'Color',Cols(id,:));
        hold on;
        end
    axis([min(minX) max(maxX) min(minY) max(maxY)].*[0.95 1.05 0.95 1.05]);

    for id=Range
        LegH(id)=Leg{id}(1);
        xV=costHistory{id,1};
        yV=costHistory{id,2};
        plot(xV(1), yV(1),'Color',Cols(id,:),'Marker','x')
        
        for n=1:sum(xV>0)-1
            p1=[xV(n) yV(n)];
            p2=[xV(n+1) yV(n+1)];
            if max(abs(diff([p1; p2])./[dX dY])) > 0.05
                arrow(p1,p2,...
                'BaseAngle',20,'Color',Cols(id,:) ,'Type','Patch');
            else
                plot(p2(1), p2(2),'Color',Cols(id,:),'Marker','o')
            end
        end
    end
    hold off;
    xlabel('Weight (kN)');    ylabel('Embodied Energy (MJ)')
%     xlabel('Weight');    ylabel('Energy')
    legend(LegH,Plots.names(Range),'Location','Best')

%     L=legend([Leg{1}(1) Leg{2}(1)],Plots.names(Range),'Location','Best')
else
    maxIt=max(cellfun(@length,costHistory(:,1)));
    its=1:maxIt;    costs=zeros(1,maxIt);   limits=0;
    for id=Range
        Ydata=costHistory{id,1};    len=length(Ydata);
        costs(1:len)=Ydata;    limits=max([ limits costs]);
        if len<maxIt
            costs(len+1:end)=Ydata(end); %fill with last values
        end
        semilogy(its,costs);    hold on
    end
   xlim([0 maxIt]); ylim([0 limits])
    xlabel('Iteration');    
    if goalID==1
        ylabel('Weight (kN)')
    else
        ylabel('Embodied Energy (MJ)')
    end

    legend(Plots.names(Range),'Location','Best')
end
grid on;    grid minor;
% fprintf('algorithm\ttime\titerations\n')
timeT=compose('%2.2fs',bestTimes(Range,1,1));

algT=Plots.names(Range);

%assuming best cost stays on the plot, which is then read from
if goalID==3
    timeT2=compose('%2.2fs',bestTimes(Range,1,2));

    algB1=cellfun(@(x) round(x(end),1),costHistory(Range,1));
    algB2=cellfun(@(x) round(x(end),1),costHistory(Range,2));
    disp(table(algT,algB1,timeT,bestTimes(Range,2,1),algB2,timeT2,bestTimes(Range,2,2),'VariableNames',...
        {'algorithm','weight','time_to_weight','iter_to_wght','energy','time_to_energy','iter_to_en'}))
else
    algB=cellfun(@(x) round(x(end),1),costHistory(Range,1));
    disp(table(algT,algB,timeT,bestTimes(Range,2,1),'VariableNames',{'algorithm','cost','time_to_best','iter_to_best'}))
end
end