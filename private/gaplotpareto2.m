function state = gaplotpareto2(options,state,flag,objectivesToPlot)
%GAPLOTPARETO Plots a Pareto front for first two objectives.
%
%   Example:
%    Create an options structure that will use GAPLOTPARETO
%    as the plot function
%     options = optimoptions('gamultiobj','PlotFcn',@gaplotpareto);

%   Copyright 2007-2016 The MathWorks, Inc.

if ~isfield(state,'Rank') || size(state.Score,2) < 2
    msg = getString(message('globaloptim:gaplotcommon:PlotFcnUnavailable','gaplotpareto'));
    title(msg,'interp','none');
    return;
end

if nargin < 4 || isempty(objectivesToPlot)
    objectivesToPlot = [1 2];
end

markers = {'og','pm','sb','dk','hr','xc'};
subPops = populationIndices(options.PopulationSize);
[~,c] = size(subPops);
for i = 1:c
    pop = subPops(:,i);
    range = pop(1):pop(2);
    tag = ['gaplotPareto',int2str(i)];
    plotFirstFront(state.Score(range,:),state.Rank(range),markers{1 + mod(i,5)},objectivesToPlot(:)',flag,tag);
    hold on
end
title(sprintf('Pareto front iteration %d',state.Generation),'interp','none')
hold off;

%---------------------------------
function plotHandle = plotFirstFront(score,nonDomRank,marker,objectivesToPlot,flag,tag)

xlabelStr = sprintf('%s %s','Objective', num2str(objectivesToPlot(1)));
ylabelStr = sprintf('%s %s','Objective', num2str(objectivesToPlot(2)));

% Check the size of score
try
    plot3d = false;
    if size(score,2) > 2 && length(objectivesToPlot) >= 3
        score = score(:,objectivesToPlot(1:3));
        plot3d = true;
    elseif size(score,2) > 1 && length(objectivesToPlot) >= 2
        score = score(:,objectivesToPlot(1:2));
    end
catch % Instead of checking indices of objectivesToPlot, set defaults and plot
    objectivesToPlot = [1;2];
    score = score(:,objectivesToPlot(1:2));
    plot3d = false;
    xlabelStr = 'Objective 1';
    ylabelStr = 'Objective 2';
    if strcmpi(flag,'init')
        warning(message('globaloptim:gaplotpareto:invalidObjectivesToPlot'));
    end
end

% Get individual from first front
minRank = 1;
xy = score((nonDomRank == minRank),:);

switch flag
    case 'init'
        if plot3d
            plotHandle = plot3(xy(:,1),xy(:,2),xy(:,3),marker);
            zlabelStr = sprintf('%s %s','Objective', num2str(objectivesToPlot(3)));
            zlabel(zlabelStr,'interp','none')
            set(gca,'View',[37 38]);
        else
            plotHandle = plot(xy(:,1),xy(:,2),marker);
        end
        set(plotHandle,'Tag',tag)
        xlabel(xlabelStr,'interp','none')
        ylabel(ylabelStr,'interp','none')
        grid on

    case {'iter', 'done'}
        plotHandle = findobj(get(gca,'Children'),'Tag',tag);
        if plot3d
            set(plotHandle,'Xdata',xy(:,1), 'Ydata',xy(:,2),'Zdata',xy(:,3));
        else
            set(plotHandle,'Xdata',score(:,1), 'Ydata',score(:,2));

%             set(plotHandle,'Xdata',xy(:,1), 'Ydata',xy(:,2));
        end
end

