function  afterOptim(type)
%evalin so that the code doesn't have to be function
evalin('base','figNum=figNum+1;');
evalin('base','times{figNum}=printTime;')

grid on;    grid minor;
%save the figure for future
evalin('base','Plots.figs{figNum} = print(''-RGBImage'');');
evalin('base',['Plots.names{figNum} =''' type ''';']);
disp(['Done ' type])

%show candidate
evalin('base','disp(''Best candidate:'')');
evalin('base','disp(x_best)');
%unwrap group to compute cost of candidate in main workspace, rounding does
%not make a difference for intGA, but called anyway to keep code shorter
evalin('base','x_best=group(round(x_best));')

evalin('base','myCost(UCs,UBs,x_best,params,10);')
evalin('base','candidate{figNum}=x_best;')
evalin('base',['constr_values=constraints(x_best,params)+1;'...
                'disp(''Constraint values:'');'...
                'disp(reshape(constr_values,8,numberElements))'])
end