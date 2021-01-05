function [ ] = BuildGroupFunction(Order)
%BuildGroupFunction Generates "group" function handle from input cell
%   Order must be a cell array, with each one of them filled with
% desired index values from the same group i.e.:
% A=[1 7 8 14];
% B=[3 5 10 12];
% C=[2 4 6];
% D=[9 11 13];
% assume x=[A B C D], then x(1) is A at specific position and so on...
% group=@(x) [x(1) x(3) x(2) x(3) x(2) x(3) x(1) x(1) x(4) x(2) x(4) x(2) x(4) x(1)];

%group function building
ngroup=length(cell2mat(Order));
aa=cell(1,ngroup);

% assume x=[A B C D], then x(1) is A at specific position and so on...
for id=1:length(Order)
    group=Order{id}; %current group indexes
    for position=1:length(group)
        %put group's number (id) in proper position deinfed in group
    	aa{group(position)}=sprintf('x(%d) ',id);
    end
end

textgroup=['group=@(x) [' cell2mat(aa) '];'];

% group=@(x) [x(1) x(3) x(2) x(3) x(2) x(3) x(1) x(1) x(4) x(2) x(4) x(2) x(4) x(1)];
evalin('base',textgroup) % use 'eval' to declare group function without modifying the code

end

