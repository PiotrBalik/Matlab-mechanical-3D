function visualize(figNum,Plots,candidate,params,continuity)
numberElements = params.numberElements;
elementNodes = params.elementNodes;
UCs=params.UCs; UBs=params.UBs;
E = params.E;   G = params.G;
nodeCoordinates=params.nodeCoordinates;

%{
if input('Display plots y/n? ','s')~='n'
    close all
    for id = 1:figNum
        figure('Name',Plots.names{id},'Position',[388 150 760 620])
        imshow(Plots.figs{id})
    end
end
%}
for num = 1:figNum
    x = candidate{num};
    
    if continuity
        error('Continuous visualization not implemented (yet)')
        for id=numberElements:-1:1
            
            SectionType=elementNodes(:,3);
            sc = SectionType (id,1);
            row=x(id);
            if sc==2 % if it's column
                
                small=find(UCs(:,27)>row,1,'last');
                upper=find(UCs(:,27)<row,1,'first');
                
                EI(id) = E* (row-UCs(small,27))*(UCs(upper,15)-UCs(small,15))/(UCs(upper,27)-UCs(small,27))+UCs(small,15);
                EA(id) = E* row;
            else
                small=find(UBs(:,27)>row,1,'last');
                upper=find(UBs(:,27)<row,1,'first');
                
                EI(id) = E* (row-UBs(small,27))*(UBs(upper,15)-UBs(small,15))/(UBs(upper,27)-UBs(small,27))+UBs(small,15);
                EA(id) = E* row;
            end
        end
    else
        for id=numberElements:-1:1
            
            SectionType=elementNodes(:,3);
            sc = SectionType (id,1);
            row=x(id);
            if sc==2 % if it's column
                EIy(id) = E* UCs(row,15); % save values to rememeber chosen row
                EIz(id) = E* UCs(row,16);
                EA(id) = E* UCs(row,27);
                GJ(id) = G* UCs(row,26); % "J is  "It" in the UCs" means J >> Torsional constant for each element
            else
                EIy(id) = E* UBs(row,15);
                EIz(id) = E* UBs(row,16);
                EA(id) = E* UBs(row,27);
                GJ(id) = G* UBs(row,26); % "J is  "It" in the UBs" means J >> Torsional constant for each element
            end
        end%for id
    end% if continuity
    
    figure('Name',['Candidate: ' Plots.names{num}])
    hold on
    forces = computeForce(EA,EIy,EIz,GJ,params);
	dM=forces.kkk;
    ds=reshape(dM,3,[]);
    tab=table(ds(1,:)',ds(2,:)',ds(3,:)','VariableNames',{'x','y','z'});
    fprintf(['\n\n' repmat('---',1,15)]);
    fprintf('\nDisplacements for candidate %d:\n',num);
    disp(tab)
% 	fprintf('%2.4G\t%2.4G\t%2.4G\n',reshape(dM,[],3));
	fprintf([repmat('---',1,15) '\n' repmat('---',1,15) '\n']);
	
	F=forces.Forces;
    fs=reshape(F,numberElements,12);
    fs2=mat2cell(fs,numberElements,ones(1,12));
    tab=table(fs2{:});
   	fprintf('Forces for candidate %d:\n',num);
    disp(tab)
% 	fprintf([repmat('%4.4G\t',[1 12])  '\n'],reshape(F,[],12));
	fprintf([repmat('-----',1,round(2.1*numberElements)) '\n\n']);
	
    magn=20; % value to increase the displacement
    for i = 1:length(elementNodes)
        
        %find shift in xy axis, from adress 1+N(i-1),
        %where i is number of node
        %N is a period (2D-two points, 3D-three points)
        Node=elementNodes(i,[1 2]); %pair of nodes
        x_shift = magn*dM(3*(Node-1)+1);
        y_shift = magn*dM(3*(Node-1)+2);
        z_shift = magn*dM(3*(Node-1)+3);
        
        %compute positions
        px=nodeCoordinates(Node,2);
        py=nodeCoordinates(Node,3);
        pz=nodeCoordinates(Node,4);
        
        %plot after displacement
        plot3(pz+z_shift,px+x_shift,py+y_shift,'bo--','MarkerSize',6)
        
        %plot before displacement
        plot3(pz,px,py,'ro--','MarkerSize',6)
    end
    grid on;    grid minor;    view([30 20]);
    xlabel('x');    ylabel('y');    zlabel('z');
    legend({'after','before'});    title(sprintf('Node shift magnified x%d times',magn))
end%for num
end%function