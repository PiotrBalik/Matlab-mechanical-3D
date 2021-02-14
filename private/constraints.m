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
E = params.E;	G = params.G;
% C1 = params.C1;
% lambdaDashLT0 = params.lambdaDashLT0;
elementNodes=params.elementNodes;
UCs = params.UCs;	UBs = params.UBs;
numberElements = params.numberElements;

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
Ned0 = forces.Ned0;	Ved0 = forces.Ved0;
MYed0 = forces.MYed0;	MZed0 = forces.MZed0;
shear = forces.shear;	moment = forces.moment;
params.ggg = forces.ggg;

for i=1:numberElements
    
    choice = x(i);
    SectionType=elementNodes(:,3);
    sc = SectionType (i,1);
    
    Ned=Ned0(i,1);		Ved=Ved0(i,1);
    MYed=MYed0(i,1);	MZed=MZed0(i,1);
    L=L0(i,1);
    shear1=shear(i,1);
%     shear2=shear(i,2);
    moment1=moment(i,1); % moment1 is M1 at the left end
%     moment2=shear1*0.25*L+moment1+Load1(i,2)*0.25*L;  % Moment at the FIRST QUARTER of the member
    moment3=shear1*0.5*L+moment1+Load1(i,2)*0.5*L*0.25*L;  % Moment at the mid of the member
%     moment4=shear1*0.75*L+moment1+Load1(i,2)*0.75*L;  % Moment at the THIRD QUARTER of the member
    moment5=moment(i,2); % moment2 is M2 at the end right

    %send values to design
    params.MYed=MYed;    params.Ned=Ned;
    params.MZed=MZed;    params.Ved=Ved;

    params.L=L;
    params.moment1=moment1;
    params.moment3=moment3;
    params.moment5=moment5;
    
    if sc==2
        UXs=UCs;
    else
        UXs=UBs;
    end
    [MR, FBRelY, FBRelZ,LTBRe1,LTBRe2,CR,SR,ED]  = parametersDesign(UXs,choice,sc,i,params);
    
    c=zeros(numberElements*7,1);
    c(1+(i-1)*8) = MR   -1;  % constraints are in form C < 0, if C < 1, then we make C - 1 < 1 - 1
    c(2+(i-1)*8) = FBRelY -1;
    c(3+(i-1)*8) = FBRelZ -1;
    c(4+(i-1)*8) = LTBRe1 -1;
    c(5+(i-1)*8) = LTBRe2 -1;
    c(6+(i-1)*8) = CR -1;
    c(7+(i-1)*8) = SR -1;
    c(8+(i-1)*8) = ED -1;
    
    ceq = [];
end
end