function [forces,displacements] = seeForces(x,params)
% global numberElements SectionType
% x - is the choice of best candidates
% so it store ind exes of UC's and UB's, to choose from

% constants reading, commented out unused ones
Load1 = params.Load1;
L0 = params.L0;

E = params.E;
G = params.G;
% C1 = params.C1;
% lambdaDashLT0 = params.lambdaDashLT0;
elementNodes=params.elementNodes;
UCs = params.UCs;
UBs = params.UBs;
numberElements = params.numberElements;
c=zeros(numberElements*7,1);
% SectionType=params.SectionType;

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
F = computeForce(EA,EIy,EIz,GJ,params);
forces=reshape(F.Forces,[],12);
displacements=reshape(F.kkk,[],3);