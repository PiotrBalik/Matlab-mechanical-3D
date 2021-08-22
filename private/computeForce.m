function NewParams = computeForce(EA,EIy,EIz,GJ,ForceParams)
% E = ForceParams.E;
% V = ForceParams.V;
% cosa = ForceParams.cosa;
% sena = ForceParams.sena;
elementDof = ForceParams.elementDof;
ll = ForceParams.ll;
loadpoints = ForceParams.loadpoints;
numberElements = ForceParams.numberElements;
GDof = ForceParams.GDof;
% elementNodes = ForceParams.elementNodes;
loadmembers = ForceParams.loadmembers;
activeDof = ForceParams.activeDof;
forceVector = ForceParams.forceVector;
memberlength = ForceParams.memberlength;
elementNodes = ForceParams.elementNodes;
xx = ForceParams.xx;
yy = ForceParams.yy;
zz = ForceParams.zz;
nodeCoordinates = ForceParams.nodeCoordinates;
DisloadDirection = ForceParams.DisloadDirection;
% steelGrade = ForceParams.steelGrade;
indexes=ForceParams.indexes;

L=zeros(1,numberElements);

stiffness = zeros(GDof);
% Tt = T;
% KgElement = zeros(12,12,numberElements);
% ki = zeros(12,12,numberElements);
% T = zeros(12,12,numberElements);
% Tt = zeros(12,12,numberElements);
%persistency improves speed
persistent a b c d KgElement TI ki T Tt
if isempty(a)
    a=zeros(3);
    b=zeros(3);
    c=zeros(3);
    d=zeros(3);
    TI = zeros(12,12);
    KgElement = zeros(12,12,numberElements);
    ki = zeros(12,12,numberElements);
    T = zeros(12,12,numberElements);
    Tt = zeros(12,12,numberElements);
end
% took constants out of loop
for i=1:numberElements

 indice=elementNodes(i,:) ; % This is only to find the length of each member below
    %   *************                                                                                                                                                                  ******************************************
    nn=length(indice);
    xa=xx(indice(2))-xx(indice(1));
    ya=yy(indice(2))-yy(indice(1));
    za=zz(indice(2))-zz(indice(1));
    length_elemen=sqrt(xa*xa+ya*ya+za*za);
    L(i)=length_elemen;
% 	LL(:,:,i)=L;
    x1=nodeCoordinates(indice(1),2);
    y1=nodeCoordinates(indice(1),3);
    z1=nodeCoordinates(indice(1),4);
    x2=nodeCoordinates(indice(2),2);
    y2=nodeCoordinates(indice(2),3);
    z2=nodeCoordinates(indice(2),4);
    
    
    k1 = EA(i)/L(i);  % 12EA/L
    k2 = 12*EIz(i)/L(i)^3; %12*EIz/L^3
    k3 = 6*EIz(i)/L(i)^2;  %6*EIz/L^2
    k4 = 4*EIz(i)/L(i);  %4*EIz/L
    k5 = 2*EIz(i)/L(i); %2*EIz/L
    k6 = 12*EIy(i)/L(i)^3; % 12*EIy/L^3
    k7 = 6*EIy(i)/L(i)^2; %6*EIy/L^2
    k8 = 4*EIy(i)/L(i); %4*EIy/L
    k9 = 2*EIy(i)/L(i); %2*EIy/L
    k10 = GJ(i)/L(i); %G*J/L
    
    
%     a=[k1 0 0; 0 k2 0; 0 0 k6];
%     b=[ 0 0 0;0 0 k3; 0 -k7 0];
%     c=[k10 0 0;0 k8 0; 0 0 k4];
%     d=[-k10 0 0;0 k9 0;0 0 k5];
    a(1,1)=k1;
    a(2,2)=k2;
    a(3,3)=k6;
%     a=[k1 0 0; 0 k2 0; 0 0 k6];
    b(2,3)=k3;
    b(3,2)=-k7;
%     b=[ 0 0 0;0 0 k3; 0 -k7 0];
    c(1,1)=k10;
    c(2,2)=k8;
    c(3,3)=k4;
%     c=[k10 0 0;0 k8 0; 0 0 k4];
    d(1,1)=-k10;
    d(2,2)=k9;
    d(3,3)=k5;
%     d=[-k10 0 0;0 k9 0;0 0 k5];
    
    %ki(:,:,i) = [a b -a b;b' c (-b)' d; (-a)' -b a -b;b' d' (-b)' c];
	KI = [a b -a b;b' c (-b)' d; (-a)' -b a -b;b' d' (-b)' c];

    
    if x1==x2 && y1==y2;
        if z2 > z1
            Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
        else
            Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
        end
    else
        CXx = (x2-x1)/L(i);
        CYx = (y2-y1)/L(i);
        CZx = (z2-z1)/L(i);
        D = sqrt(CXx*CXx + CYx*CYx);
        CXy = -CYx/D;
        CYy = CXx/D;
        CZy = 0;
        CXz = -CXx*CZx/D;
        CYz = -CYx*CZx/D;
        CZz = D;
        Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
    end
    
    
%     T(:,:,i) = [Lambda zeros(3,9); zeros(3) Lambda zeros(3,6);
%         zeros(3,6) Lambda zeros(3);zeros(3,9) Lambda];
% Around 3x faster with indexed Lambda
TI(1:3,1:3) = Lambda;
TI(4:6,4:6) = Lambda;
TI(7:9,7:9) = Lambda;
TI(10:12,10:12) = Lambda;

T(:,:,i)=TI;
Tt(:,:,i)=TI';
%Tt(:,:,i)=T(:,:,i)';
KgElement(:,:,i)=TI'*KI*TI;
%KgElement(:,:,i)=Tt(:,:,i)*KI*T(:,:,i);
ki(:,:,i) = KI;
stiffness(elementDof(i,:),elementDof(i,:))= stiffness(elementDof(i,:),elementDof(i,:))+KgElement(:,:,i);
%     stiffness(elementDof(i,:),elementDof(i,:))= stiffness(elementDof(i,:),elementDof(i,:))+Tt(:,:,i)*ki(:,:,i)*T(:,:,i);

end

% Displacement Calculations

%!!!------line below moved to precomputation----!!!
% activeDof=setdiff([1:GDof]',[prescribedDof]);
% copy of working functions - requires sorted data and proper classes

U=stiffness(activeDof,activeDof)\forceVector(activeDof);
displacements=zeros(GDof,1);
displacements(activeDof)=U;
ggg=displacements(1:6:end);

%workaround to remove every 4th,5th,6th element: by indexation (fastest)
kkk=displacements(indexes);
% kkk=displacements(mod(1:length(displacements),4)~=0);
% kkk=displacements(setdiff(1:length(displacements),3:3:length(displacements)));

% Internal forces of each element

distributedloadReactions=zeros(numberElements,12);
ss=size(loadmembers);

for i=1:ss(1)
    
    rr=loadmembers(i,1);   % to find the element number
    
    L=memberlength(rr,1);
    
    P=loadmembers(loadmembers(:,1)==rr,2);
    
    % Load Direction
    
    % Xy ... member X and load in y direction
    
        if DisloadDirection(i,1)==1  %condition for the member direction and load direction,
    
    distributedloadReactions(rr,2)=-P*L/2; % IT WAS -
    distributedloadReactions(rr,8)=-P*L/2; % IT WAS -
    
    distributedloadReactions(rr,6)=-P*L^2/12; % CHANGED TRIAL ---***----****-----**** % IT WAS -
    distributedloadReactions(rr,12)=P*L^2/12; % CHANGED TRIAL ---***----****-----**** % IT WAS +
        end
    
    % Xz ... member X and load in z direction
    
        if DisloadDirection(i,1)==2  %condition for the member direction and load direction,

    distributedloadReactions(rr,3)=-P*L/2;
    distributedloadReactions(rr,9)=-P*L/2;
    
    distributedloadReactions(rr,5)=-P*L^2/12; % IT WAS +
    distributedloadReactions(rr,11)=P*L^2/12; % IT WAS -

        end
        
    % Yx ... member Y and load in x direction
    
        if DisloadDirection(i,1)==3  %condition for the member direction and load direction,

    distributedloadReactions(rr,1)=P*L/2;  % changed , which was 1
    distributedloadReactions(rr,7)=P*L/2;  % changed , which was 7
    
    distributedloadReactions(rr,6)=P*L^2/12; % CHANGED TRIAL ---***----****-----**** % IT WAS +
    distributedloadReactions(rr,12)=-P*L^2/12; % CHANGED TRIAL ---***----****-----**** % IT WAS -
        end
    
    % Yz ... member Y and load in z direction
    
        if DisloadDirection(i,1)==4  %condition for the member direction and load direction,

    distributedloadReactions(rr,3)=-P*L/2; % IT WAS -
    distributedloadReactions(rr,9)=-P*L/2; % IT WAS -
    
    distributedloadReactions(rr,4)=-P*L^2/12;
    distributedloadReactions(rr,10)=P*L^2/12;
    
        end
        
    % Zx ... member Z and load in x direction
    
        if DisloadDirection(i,1)==5  %condition for the member direction and load direction,

    distributedloadReactions(rr,2)=-P*L/2; % changed , which was 1
    distributedloadReactions(rr,8)=-P*L/2; % changed , which was 7
    
    distributedloadReactions(rr,5)=-P*L^2/12;
    distributedloadReactions(rr,11)=P*L^2/12;
        end
        
    % Zy ... member Z and load in y direction

        if DisloadDirection(i,1)==6  %condition for the member direction and load direction,

    distributedloadReactions(rr,2)=-P*L/2; % IT WAS - >>>2 -
    distributedloadReactions(rr,8)=-P*L/2; % IT WAS - >>>8 -
    
    distributedloadReactions(rr,6)=-P*L^2/12;   %<<<<<<<<< IT WAS 4 >>>6 -
    distributedloadReactions(rr,12)=P*L^2/12; %<<<<<<<<< IT WAS 10 >>>>12 +
        end
    
    
end
%Distributed load on each member
%Point load on each member

pointloadReactions=zeros(numberElements,12);
jj=size(loadpoints);

for i=1:jj(1)
    
 L=memberlength(i,1);
    
    rr=loadpoints(i,1); % to find the elemnet number from matrix 
    
    L=memberlength(rr,1);
    
    P=loadpoints(loadpoints(:,1)==rr,2);
    
    a=loadpoints(loadpoints(:,1)==rr,3);
    
    b=L-a;
    
    % *** Load direction conditions *****************
    
    
    % Xy ... member X and load in y direction
    
    
    if PointloadDirection(i,1)==1;  %        %condition for the member direction and load direction,
        
        
        if a>b || a<b
            
            pointloadReactions(rr,2)=-(P*b^2/L^3)*(3*a+b);
            pointloadReactions(rr,8)=-(P*a^2/L^3)*(a+3*b);
            
            pointloadReactions(rr,6)=-(P*a*b^2/L^2); % CHANGED TRIAL ---***----****-----****
            pointloadReactions(rr,12)=(P*a^2*b/L^2); % CHANGED TRIAL ---***----****-----****
        end
        
        if a==b;
            
            pointloadReactions(rr,2)=-P/2;
            pointloadReactions(rr,8)=-P/2;
            
            pointloadReactions(rr,6)=-P*L/8; % CHANGED TRIAL ---***----****-----****
            pointloadReactions(rr,12)=P*L/8; % CHANGED TRIAL ---***----****-----****
        end
        
    end
    
    
    % Xz ... member X and load in z direction
    
    if PointloadDirection(i,1)==2;  %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointloadReactions(rr,3)=-(P*b^2/L^3)*(3*a+b);
            pointloadReactions(rr,9)=-(P*a^2/L^3)*(a+3*b);
            
            pointloadReactions(rr,5)=(P*a*b^2/L^2);
            pointloadReactions(rr,11)=-(P*a^2*b/L^2);
        end
        
        if a==b;
            
            pointloadReactions(rr,3)=-P/2;
            pointloadReactions(rr,9)=-P/2;
            
            pointloadReactions(rr,5)=P*L/8;
            pointloadReactions(rr,11)=-P*L/8;
        end
        
    end
    
    
    % Yx ... member Y and load in x direction
    
    if PointloadDirection(i,1)==3;  %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointloadReactions(rr,2)=-(P*b^2/L^3)*(3*a+b);  % changed , which was 1
            pointloadReactions(rr,8)=-(P*a^2/L^3)*(a+3*b);  % changed , which was 7
            
            pointloadReactions(rr,6)=(P*a*b^2/L^2); % CHANGED TRIAL ---***----****-----****
            pointloadReactions(rr,12)=-(P*a^2*b/L^2); % CHANGED TRIAL ---***----****-----****
        end
        
        if a==b;
            
            pointloadReactions(rr,2)=-P/2;
            pointloadReactions(rr,8)=-P/2;
            
            pointloadReactions(rr,6)=P*L/8; % CHANGED TRIAL ---***----****-----****
            pointloadReactions(rr,12)=-P*L/8; % CHANGED TRIAL ---***----****-----****
        end
        
        
    end
    
    
    % Yz ... member Y and load in z direction
    
    if PointloadDirection(i,1)==4;  %condition for the member direction and load direction,
        
        
        if a>b || a<b
            
            pointloadReactions(rr,3)=-(P*b^2/L^3)*(3*a+b);
            pointloadReactions(rr,9)=-(P*a^2/L^3)*(a+3*b);
            
            pointloadReactions(rr,4)=-(P*a*b^2/L^2);
            pointloadReactions(rr,10)=(P*a^2*b/L^2);
        end
        
        if a==b;
            
            pointloadReactions(rr,3)=-P/2;
            pointloadReactions(rr,9)=-P/2;
            
            pointloadReactions(rr,4)=-P*L/8;
            pointloadReactions(rr,10)=P*L/8;
        end
        
        
    end
    
    % Zx ... member Z and load in x direction
    
    if PointloadDirection(i,1)==5;  %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointloadReactions(rr,2)=-(P*b^2/L^3)*(3*a+b); % changed , which was 2
            pointloadReactions(rr,8)=-(P*a^2/L^3)*(a+3*b); % changed , which was 8
            
            pointloadReactions(rr,5)=-(P*a*b^2/L^2);
            pointloadReactions(rr,11)=(P*a^2*b/L^2);
        end
        
        if a==b;
            
            pointloadReactions(rr,2)=-P/2;
            pointloadReactions(rr,8)=-P/2;
            
            pointloadReactions(rr,5)=-P*L/8;
            pointloadReactions(rr,11)=P*L/8;
        end
    end
    
    % Zy ... member Z and load in y direction
    
    if PointloadDirection(i,1)==6 ; %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointloadReactions(rr,2)=-(P*b^2/L^3)*(3*a+b);
            pointloadReactions(rr,8)=-(P*a^2/L^3)*(a+3*b);
            
            pointloadReactions(rr,4)=(P*a*b^2/L^2);
            pointloadReactions(rr,10)=-(P*a^2*b/L^2);
        end
        
        if a==b;
            
            pointloadReactions(rr,2)=-P/2;
            pointloadReactions(rr,8)=-P/2;
            
            pointloadReactions(rr,4)=P*L/8; % 4
            pointloadReactions(rr,10)=-P*L/8; % 10
        end
        
    end
    
    
end

ww=size(elementDof);   % to find the ww1 which is the total number of Dof
ww1=ww(1)*ww(2); %ww is the size of all Dof of all elements,(all elemens=6*number of elements)
indxx=reshape(elementDof,[ww1 1]);

mm=reshape(distributedloadReactions,[ww1 1]);
nn=reshape(pointloadReactions,[ww1 1]);
indxx(:,2)=mm;
indxx(:,3)=nn;
                                % HERE I STOPPPPPPPPED 
                                
                                
for i=12:-1:1 % -1:1 is implicit preallocation from reverse indexing, right, it was not there in Frame3D File :)
    Reactions(i,1)=sum(indxx(indxx(:,1)==i,2));
    Reactions(i,2)=sum(indxx(indxx(:,1)==i,3));
end
% small cleanup:
% ElementReactions=zeros(12,2);
% ElementReactions(:,3)=Reactions(:,1)+Reactions(:,2);

% ElementReactions=Reactions(:,1)+Reactions(:,2); unused?
% calculations loop
Reactions11=pointloadReactions+distributedloadReactions;
Reaction22=Reactions11';

% ElementReactions(:,:,i)=Reaction22(:,i);
DofEachElement = zeros(ww(2),1,numberElements);
Di=DofEachElement;
internalForces=Di;

for i=1:numberElements
    
     ElementReactions=Reaction22(:,i);
    
    DofEachElement(:,:,i)=elementDof(i,:)';
    
    Di(:,:,i)=displacements(DofEachElement(:,:,i),1);
    
    internalForces(:,:,i)=ki(:,:,i)*T(:,:,i)*Di(:,:,i)+ElementReactions;
end

% to make memberConnectivities in one column then transfer them to excel sheet
% Q=numberElements*2;  % TO FIND THE TOTAL NUMBER OF NODES FOR ALL MEMBERS, ( total number of members * 2 because each member has 2 nodes)
% memberConnectivities=reshape((elementNodes(:,1:2))',[Q 1]);

% to make internal forces in one column then trensfer to excel
ww=size(elementDof);   % to find the ww1 which is the total number of Dof
ww1=ww(1)*ww(2);

memberforces=reshape(internalForces,[ww1 1]);

% epsilon=sqrt(235/steelGrade);

% G = E/(2*(1+V));% shear modulus

% lambda1=93.9*epsilon;

% Reading the AXIAL (Ned), MOMENT (Med),SHEAR(Ved)----------------------

AxialNed=memberforces(1:6:end);
% AxialNed= xlsread('Input3D.xlsm','output', 'O5:O100');  % READING THE AXIAL FORCE (Ned)AND TAKING THE HIGEAST FOR EACH MEMBER
Ned1=abs(AxialNed(1:2:end));
Ned1(:,2)=abs(AxialNed(2:2:end));
Ned0=(max(Ned1'))';

ShearVed=memberforces(2:6:end);
% ShearVed= xlsread('Input3D.xlsm','output', 'P5:P100');  % READING THE SHEAR FORCE (Ved)AND TAKING THE HIGEAST FOR EACH MEMBER
Ved1=abs(ShearVed(1:2:end));
Ved1(:,2)=abs(ShearVed(2:2:end));
Ved0=(max(Ved1'))';

FZ=memberforces(3:6:end);
% MomentMed= xlsread('Input3D.xlsm','output', 'Q5:Q100');  % READING THE Bending moment (Med) AND TAKING THE HIGEAST FOR EACH MEMBER
FZed1=abs(FZ(1:2:end));
FZed1(:,2)=abs(FZ(2:2:end));
FZed0=(max(FZed1'))';
    

MX=memberforces(4:6:end);
% MomentMed= xlsread('Input3D.xlsm','output', 'Q5:Q100');  % READING THE Bending moment (Med) AND TAKING THE HIGEAST FOR EACH MEMBER
MXed1=abs(MX(1:2:end));
MXed1(:,2)=abs(MX(2:2:end));
MXedx0=(max(MXed1'))';


MY=memberforces(5:6:end);
% MomentMed= xlsread('Input3D.xlsm','output', 'Q5:Q100');  % READING THE Bending moment (Med) AND TAKING THE HIGEAST FOR EACH MEMBER
MYd1=abs(MY(1:2:end));
MYed1(:,2)=abs(MY(2:2:end));
MYed0=(max(MYed1'))';

MZ=memberforces(6:6:end);
% MomentMed= xlsread('Input3D.xlsm','output', 'Q5:Q100');  % READING THE Bending moment (Med) AND TAKING THE HIGEAST FOR EACH MEMBER
MZed1=abs(MZ(1:2:end));
MZed1(:,2)=abs(MZ(2:2:end));
MZed0=(max(MZed1'))';

%%% to find matrix for the load on each member (match and fill in one matrix""named Load1"") %%%

shear=ShearVed(1:2:end);
shear(:,2)=ShearVed(2:2:end);
moment=MZ(1:2:end);
moment(:,2)=MZ(2:2:end);


    NewParams.Forces=memberforces;
    NewParams.Ned0 = Ned0;
    NewParams.Ved0 = Ved0;
    NewParams.FZed0 = FZed0;
    NewParams.MXedx0 = MXedx0;
    NewParams.MYed0 = MYed0;
    NewParams.MZed0 = MZed0;

    NewParams.shear = shear;
    NewParams.moment = moment;
    NewParams.ggg=ggg;
    NewParams.kkk=kkk;

end