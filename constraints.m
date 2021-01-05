function [c,ceq] = constraints(x,params)
% global numberElements SectionType
% x - is the choice of best candidates
% so it store ind exes of UC's and UB's, to choose from

% constants reading, commented out unused ones
Load1 = params.Load1;
L0 = params.L0;
% gammaM0 = params.gammaM0;
% gammaM1 = params.gammaM1;
% steelGrade = params.steelGrade;
% epsilon = params.epsilon;
% ita = params.ita;
% beta = params.beta;
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
forces = computeForce(EA,EIy,EIz,GJ,params);
Ned0 = forces.Ned0;
Ved0 = forces.Ved0;
MYed0 = forces.MYed0;
MZed0 = forces.MZed0;
shear = forces.shear;
moment = forces.moment;
params.ggg = forces.ggg;

for i=1:numberElements
    
    choice = x(i);
    SectionType=elementNodes(:,3);
    sc = SectionType (i,1);
    
    Ned=Ned0(i,1);
    Ved=Ved0(i,1);
    MYed=MYed0(i,1);
    MZed=MZed0(i,1);
    L=L0(i,1);
    shear1=shear(i,1);
%     shear2=shear(i,2);
    moment1=moment(i,1); % moment1 is M1 at the left end
%     moment2=shear1*0.25*L+moment1+Load1(i,2)*0.25*L;  % Moment at the FIRST QUARTER of the member
    moment3=shear1*0.5*L+moment1+Load1(i,2)*0.5*L*0.25*L;  % Moment at the mid of the member
%     moment4=shear1*0.75*L+moment1+Load1(i,2)*0.75*L;  % Moment at the THIRD QUARTER of the member
    moment5=moment(i,2); % moment2 is M2 at the end right

    %send values to design
    params.Ned=Ned;
    params.Ved=Ved;
    params.MYed=MYed;
    params.MZed=MZed;

    params.L=L;
    params.moment1=moment1;
    params.moment3=moment3;
    params.moment5=moment5;
    
    if sc==2
        [MRc, FBRColumnY, FBRColumnZ,LTBRC1,LTBRC2,Crc,SRc,VD]  = parametersDesign(UCs,choice,sc,i,params);
        c(1+(i-1)*8)= MRc  -1;  % constraints are in form C < 0, if C < 1, then we make C - 1 < 1 - 1
        c(2+(i-1)*8)= FBRColumnY  -1;
         c(3+(i-1)*8)=  FBRColumnZ -1;
        c(4+(i-1)*8)= LTBRC1  -1;
        c(5+(i-1)*8)= LTBRC2  -1;
        c(6+(i-1)*8)= Crc  -1;
        c(7+(i-1)*8)=  SRc -1;
        c(8+(i-1)*8)= VD - 1;
    else
        [MRB, FBRBeamY, FBRBeamZ,LTBRB1,LTBRB2,CrB,SRB,HD]  = parametersDesign(UBs,choice,sc,i,params);
        c(1+(i-1)*8)= MRB  -1;
        c(2+(i-1)*8)=  FBRBeamY -1;
        c(3+(i-1)*8)=  FBRBeamZ -1;
        c(4+(i-1)*8)=  LTBRB1 -1;
        c(5+(i-1)*8)=  LTBRB2 -1;
        c(6+(i-1)*8)=  CrB -1;
        c(7+(i-1)*8)=  SRB -1;
        c(8+(i-1)*8)= HD - 1;

    end
    
    ceq = [];
end
end