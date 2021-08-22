function  afterOptim(afterParams)

is2D=afterParams.is2D;
type=afterParams.type;
multiObj=afterParams.multiObj;
grouping=afterParams.grouping;
continuity=afterParams.continuity;
Nobjectives=afterParams.Nobjectives;

%evalin so that the can access base workspace
%and main doesn't need to be a function
grid on;    grid minor;
evalin('base','figNum=figNum+1;');
evalin('base','times(figNum)=printTime;')

if multiObj==1
    for nObj=1:Nobjectives % now we can have 100 objectives at once :)
    evalin('base',sprintf('bestTimes(figNum,1,%d)=history.bestTime(%d);',nObj,nObj));
    evalin('base',sprintf('bestTimes(figNum,2,%d)=history.bestIt(%d);',nObj,nObj));
    evalin('base',sprintf('costHistory{figNum,%d}=history.Costs(:,%d);',nObj,nObj));
    end
else
    evalin('base','Line=get(gca,''Children'');Line=Line(end);');
    evalin('base','[~,p]=min(Line.YData);');
    evalin('base','pmax=Line.XData(end);');
    evalin('base','bestTimes(figNum,1,1)=times(figNum)*p/pmax;');
    evalin('base','bestTimes(figNum,2,1)=p;');
    evalin('base','costHistory{figNum,1}=Line.YData;');
end
%save the figure for future
evalin('base','Plots.figs{figNum} = print(''-RGBImage'');');
evalin('base',['Plots.names{figNum} =''' type ''';']);
disp(['Done ' type])

if continuity==0 %discrete
    evalin('base','x_best=round(x_best);')
end

%show candidate
evalin('base','disp(''Best candidate:'')');
evalin('base','disp(x_best)');
evalin('base','candidate{figNum}=x_best;')
%unwrap group to compute cost of candidate in main workspace
if grouping
	evalin('base','x_best=group(x_best);');
end

if continuity
    evalin('base','myCost(UCCs,UBBs,x_best,params,10);')
else
    evalin('base','myCost(UCs,UBs,x_best,params,10);')
end

if is2D
	evalin('base',['constr_values=constraints(x_best,params)+1;'...
		'disp(''Constraint values:'');'...
		'disp(reshape(constr_values,7,numberElements))'])
else
	evalin('base',['constr_values=constraints(x_best,params)+1;'...
		'disp(''Constraint values:'');'...
		'disp(reshape(constr_values,8,numberElements))'])
end