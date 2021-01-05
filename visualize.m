function visualize(figNum,Plots,candidate,params,continuity)
numberElements = params.numberElements;
elementNodes = params.elementNodes;
UCs=params.UCs;
UBs=params.UBs;
E = params.E;
G = params.G;
nodeCoordinates=params.nodeCoordinates;

if input('Visualize plots y/n? ','s')~='n'
    close all
    for id = 1:figNum
        figure('Name',Plots.names{id},'Position',[388 150 760 620])
        imshow(Plots.figs{id})
    end
end


for num = 1:figNum
    x = candidate{num};
    
    if continuity
        error('Not implemented (yet)')
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
                EIz(id) = E* UCs(row,16); % save values to rememeber chosen row
                EA(id) = E* UCs(row,27);
                GJ(id) = G* UCs(row,26); % "J is  "It" in the UCs" means J >> Torsional constant for each element
                
            else
                EIy(id) = E* UBs(row,15);
                EIz(id) = E* UBs(row,16); % save values to rememeber chosen row
                EA(id) = E* UBs(row,27);
                GJ(id) = G* UBs(row,26); % "J is  "It" in the UBs" means J >> Torsional constant for each element
            end
        end
    end
    
    figure('Name',['Candidate: ' Plots.names{num}])
    hold on
    forces = computeForce(EA,EIy,EIz,GJ,params);
    magnitude=20; % value to increase the displacement
    for i = 1:length(elementNodes)
        
        %find shift in xy axis, from adress 1+N(i-1),
        %where i is number of node
        %N is a period (2D-two points, 3D-three points)
        x_shift = magnitude*forces.kkk(1+3*(elementNodes(i,[1 2])-1));
        y_shift = magnitude*forces.kkk(2+3*(elementNodes(i,[1 2])-1));
        z_shift = magnitude*forces.kkk(3*elementNodes(i,3));
%         z_shift=0;
        
        %compute positions
        px=nodeCoordinates(elementNodes(i,[1 2]),2);
        py=nodeCoordinates(elementNodes(i,[1 2]),3);
        pz=nodeCoordinates(elementNodes(i,[1 2]),4);
        
        %plot after displacement
        plot3(pz+z_shift,px+x_shift,py+y_shift,'bo--','MarkerSize',6)
        
        %plot before displacement
        plot3(pz,px,py,'ro--','MarkerSize',6)
    end
    grid on
    grid minor
    view([30 20])
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend({'after','before'})
    title(sprintf('Node shift magnified by x%d',magnitude))
    % for i = 1:length(nodeCoordinates)
    %     text_add = sprintf('x = %2.3f mm\n y = %2.3f mm',10*forces.kkk((1+2*(i-1):2*i)));
    %     text(nodeCoordinates(i,2),nodeCoordinates(i,3)-0.2,text_add)
    % end
end
end