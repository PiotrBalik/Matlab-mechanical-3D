wclear all

nodeCoordinates= xlsread('Input3D.xlsm', 'A21:D500');
elementNodes=xlsread('Input3D.xlsm', 'F21:H500');
supportNode=xlsread('Input3D.xlsm', 'R21:R500'); % reading the supports number
loadjoints=xlsread('Input3D.xlsm', 'W21:Z500');   % Nodal load
loadmembers=xlsread('Input3D.xlsm', 'AE21:AF500');  % Distributed load
loadpoints=xlsread('Input3D.xlsm', 'AK21:AM500');  % point load
numberElements=xlsread('Input3D.xlsm', 'E6:E6');
DisloadDirection=xlsread('Input3D.xlsm', 'AD21:AD500');
PointloadDirection=xlsread('Input3D.xlsm', 'AJ21:AJ500');

steelGrade= xlsread('Input3D.xlsm', 'J5:J5');
gammaM0=xlsread('Input3D.xlsm', 'J7:J7'); % Partial factor (resistance of cross-section whatever the class is)
gammaM1=xlsread('Input3D.xlsm', 'J8:J8'); % Partial factor (resistance of members to instability)
V= xlsread('Input3D.xlsm', 'J9:J9');    % READING THE Passion's ratio in elastic stage (V)
E= xlsread('Input3D.xlsm', 'N21:N21');    % READING THE MODULUS OF ELASTICITY (E)
beta= xlsread('Input3D.xlsm', 'J10:J10');
ita= xlsread('Input3D.xlsm', 'J6:J6');
lambdaDashLT0=xlsread('Input3D.xlsm', 'J11:J11');
C1=xlsread('Input3D.xlsm', 'J12:J12'); %factor debends on shape of Bending moment

pointloadADirection=xlsread('Input3D.xlsm', 'AE21:AE500');
LiveLoad=xlsread('Input3D.xlsm', 'AG21:AG500'); % NEEDED ONLY FOR VERTICAL DEFLECTION************************************++++

elementDof=[];
elementDof(:,1)=elementNodes(:,1)*6-5;
elementDof(:,2)=elementNodes(:,1)*6-4;
elementDof(:,3)=elementNodes(:,1)*6-3;
elementDof(:,4)=elementNodes(:,1)*6-2;
elementDof(:,5)=elementNodes(:,1)*6-1;
elementDof(:,6)=elementNodes(:,1)*6;
elementDof(:,7)=elementNodes(:,2)*6-5;
elementDof(:,8)=elementNodes(:,2)*6-4;
elementDof(:,9)=elementNodes(:,2)*6-3;
elementDof(:,10)=elementNodes(:,2)*6-2;  
elementDof(:,11)=elementNodes(:,2)*6-1;
elementDof(:,12)=elementNodes(:,2)*6;

%%% To read PRESCRIBED DOF (the supports degree of freedom to find the displacement and reactions later) %%%%%%%%

prescribedDof0=[];
prescribedDof0(:,1)=supportNode(:,1)*6-5;
prescribedDof0(:,2)=supportNode(:,1)*6-4;
prescribedDof0(:,3)=supportNode(:,1)*6-3;
prescribedDof0(:,4)=supportNode(:,1)*6-2;
prescribedDof0(:,5)=supportNode(:,1)*6-1;
prescribedDof0(:,6)=supportNode(:,1)*6;


prescribedDof00=prescribedDof0';

SN=size(supportNode);   % SN is the number of notes that have supports
SD=SN(1)*6; % Total degree of freedom at the nodes which have supports / 6 is the Dof each node
prescribedDof=reshape(prescribedDof00,[SD 1]);

numberNodes=size(nodeCoordinates,1);
GDof=6*numberNodes;

xx=nodeCoordinates(:,2);
yy=nodeCoordinates(:,3);
zz=nodeCoordinates(:,4);

%Length and angle of each member
SectionType=elementNodes(:,3);

memberlength=zeros(numberElements,3);

for i=1:numberElements
    
    indice=elementNodes(i,:) ; % This is only to find the length of each member below
    %   *************                                                                                                                                                                  ******************************************
    nn=length(indice);
    xa=xx(indice(2))-xx(indice(1));
    ya=yy(indice(2))-yy(indice(1));
    za=zz(indice(2))-zz(indice(1));
    length_elemen=sqrt(xa*xa+ya*ya+za*za);
    ll=length_elemen;
    memberlength(i,1)=ll;
    
    
end

distributedload=zeros(numberElements,12);
ss=size(loadmembers);

for i=1:ss(1)
    
    rr=loadmembers(i,1);   % to find the element number
    
    L=memberlength(rr,1);
    
    P=loadmembers(loadmembers(:,1)==rr,2);
    
    % Load Driction
    
    % Xy ... member X and load in y direction
    
        if DisloadDirection(i,1)==1  %condition for the member direction and load direction,
    
    distributedload(rr,2)=-P*L/2;
    distributedload(rr,8)=-P*L/2;
    
    distributedload(rr,6)=-P*L^2/12;  % I changed sign 
    distributedload(rr,12)=P*L^2/12; % I changed sign
        end
    
    % Xz ... member X and load in z direction
    
        if DisloadDirection(i,1)==2  %condition for the member direction and load direction,

    distributedload(rr,3)=-P*L/2;
    distributedload(rr,9)=-P*L/2;
    
    distributedload(rr,5)=P*L^2/12;
    distributedload(rr,11)=-P*L^2/12;

        end
        
    % Yx ... member Y and load in x direction
    
        if DisloadDirection(i,1)==3  %condition for the member direction and load direction,

    distributedload(rr,1)=-P*L/2; 
    distributedload(rr,7)=-P*L/2;
    
    distributedload(rr,6)=P*L^2/12;  % Changed this sign temporary *****-***---*****-***
    distributedload(rr,12)=-P*L^2/12;  % Changed this sign temporary ***---*****-***-*****
        end
    
    % Yz ... member Y and load in z direction
    
        if DisloadDirection(i,1)==4  %condition for the member direction and load direction,

    distributedload(rr,3)=-P*L/2;
    distributedload(rr,9)=-P*L/2;
    
    distributedload(rr,4)=-P*L^2/12;
    distributedload(rr,10)=P*L^2/12;
    
        end
        
    % Zx ... member Z and load in x direction
    
        if DisloadDirection(i,1)==5  %condition for the member direction and load direction,

    distributedload(rr,1)=-P*L/2;
    distributedload(rr,7)=-P*L/2;
    
    distributedload(rr,5)=-P*L^2/12;
    distributedload(rr,11)=P*L^2/12;
        end
        
    % Zy ... member Z and load in y direction

        if DisloadDirection(i,1)==6  %condition for the member direction and load direction,

    distributedload(rr,2)=-P*L/2;
    distributedload(rr,8)=-P*L/2;
    
    distributedload(rr,4)=P*L^2/12;  % it was +  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    distributedload(rr,10)=-P*L^2/12;  % it was -  
        end
    
    
        
    
    
end

% distributedloadReaction=distributedload;
distributedloadAction=(distributedload)*-1;  %%% change directions in order to get the actual load on structure


%Nodal load on each joint

Nodalload=zeros(numberNodes,6);
% L=5;
vv=size(loadjoints);

for i=1:vv(1)
    
    ww=loadjoints(i,1);
    N1=loadjoints(loadjoints(:,1)==ww,2);
    N2=loadjoints(loadjoints(:,1)==ww,3);
    N3=loadjoints(loadjoints(:,1)==ww,4);
    
    Nodalload(ww,1)=N1;
    Nodalload(ww,2)=N2;
    Nodalload(ww,3)=N3;
    
end


pointload=zeros(numberElements,12);
% L=5;
jj=size(loadpoints);

for i=1:jj(1)
    
%     L=memberlength(i,1);
    
    rr=loadpoints(i,1); % to find the elemnet number from matrix 
    
    L=memberlength(rr,1);
    
    P=loadpoints(loadpoints(:,1)==rr,2);
    
    a=loadpoints(loadpoints(:,1)==rr,3);
    
    
    b=L-a;
    
    % *** Load direction conditions *****************
    
    
    % Xy ... member X and load in y direction
    
    
    if PointloadDirection(i,1)==1;  %        %condition for the member direction and load direction,
        
        
        
        if a>b || a<b
            
            pointload(rr,2)=-(P*b^2/L^3)*(3*a+b);
            pointload(rr,8)=-(P*a^2/L^3)*(a+3*b);
            
            pointload(rr,6)=-(P*a*b^2/L^2); % I HAVE CHANGED THIS SIGN
            pointload(rr,12)=(P*a^2*b/L^2); % I HAVE CHANGED THIS SIGN
            
        end
        
        if a==b;
            
            pointload(rr,2)=-P/2;
            pointload(rr,8)=-P/2;
            
            pointload(rr,6)=-P*L/8; % I HAVE CHANGED THIS SIGN
            pointload(rr,12)=P*L/8; % I HAVE CHANGED THIS SIGN
        end
        
    end
    
    
    % Xz ... member X and load in z direction
    
    if PointloadDirection(i,1)==2;  %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointload(rr,3)=-(P*b^2/L^3)*(3*a+b);
            pointload(rr,9)=-(P*a^2/L^3)*(a+3*b);
            
            pointload(rr,5)=(P*a*b^2/L^2);
            pointload(rr,11)=-(P*a^2*b/L^2);
            
        end
        
        if a==b;
            
            pointload(rr,3)=-P/2;
            pointload(rr,9)=-P/2;
            
            pointload(rr,5)=P*L/8;
            pointload(rr,11)=-P*L/8;
        end
        
    end
    
    
    % Yx ... member Y and load in x direction
    
    if PointloadDirection(i,1)==3;  %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointload(rr,1)=-(P*b^2/L^3)*(3*a+b);
            pointload(rr,7)=-(P*a^2/L^3)*(a+3*b);
            
            pointload(rr,6)=(P*a*b^2/L^2); % I have changed this sign
            pointload(rr,12)=-(P*a^2*b/L^2); % I have changed this sign
            
        end
        
        if a==b;
            
            pointload(rr,1)=-P/2;
            pointload(rr,7)=-P/2;
            
            pointload(rr,6)=P*L/8; % I have changed this sign
            pointload(rr,12)=-P*L/8; % I have changed this sign
        end
        
        
    end
    
    
    % Yz ... member Y and load in z direction
    
    if PointloadDirection(i,1)==4;  %condition for the member direction and load direction,
        
        
        if a>b || a<b
            
            pointload(rr,3)=-(P*b^2/L^3)*(3*a+b);
            pointload(rr,9)=-(P*a^2/L^3)*(a+3*b);
            
            pointload(rr,4)=-(P*a*b^2/L^2);
            pointload(rr,10)=(P*a^2*b/L^2);
            
        end
        
        if a==b;
            
            pointload(rr,3)=-P/2;
            pointload(rr,9)=-P/2;
            
            pointload(rr,4)=-P*L/8;
            pointload(rr,10)=P*L/8;
        end
        
        
    end
    
    % Zx ... member Z and load in x direction
    
    if PointloadDirection(i,1)==5;  %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointload(rr,1)=-(P*b^2/L^3)*(3*a+b);
            pointload(rr,7)=-(P*a^2/L^3)*(a+3*b);
            
            pointload(rr,5)=-(P*a*b^2/L^2);
            pointload(rr,11)=(P*a^2*b/L^2);
            
        end
        
        if a==b;
            
            pointload(rr,1)=-P/2;
            pointload(rr,7)=-P/2;
            
            pointload(rr,5)=-P*L/8;
            pointload(rr,11)=P*L/8;
        end
    end
    
    % Zy ... member Z and load in y direction
    
    if PointloadDirection(i,1)==6 ; %condition for the member direction and load direction,
        
        if a>b || a<b
            
            pointload(rr,2)=-(P*b^2/L^3)*(3*a+b);
            pointload(rr,8)=-(P*a^2/L^3)*(a+3*b);
            
            pointload(rr,4)=(P*a*b^2/L^2);
            pointload(rr,10)=-(P*a^2*b/L^2);
            
        end
        
        if a==b;
            
            pointload(rr,2)=-P/2;
            pointload(rr,8)=-P/2;
            
            pointload(rr,4)=P*L/8;
            pointload(rr,10)=-P*L/8;
        end
        
    end
    
    
    
    
end
    

% pointloadReaction=pointload;
pointloadAction=(pointload);  %%% change directions in order to get the actual load on structure

ll=memberlength(:,1);

forceVector=zeros(GDof,1);


ww=size(elementDof);   % to find the ww1 which is the total number of Dof
ww1=ww(1)*ww(2); %ww is the size of all Dof of all elements,(all elemens=6*number of elements)

indx=reshape(elementDof,[ww1 1]);

aa=reshape(distributedloadAction,[ww1 1]);
bb=reshape(pointloadAction,[ww1 1]);
indx(:,2)=aa;
indx(:,3)=bb;

for i=1:GDof
    cc(i,1)=sum(indx(indx(:,1)==i,2));
    cc(i,2)=sum(indx(indx(:,1)==i,3));
    
end

cc(:,3)=reshape(Nodalload',[GDof 1]);
cc(:,4)=cc(:,1)+cc(:,2)+cc(:,3);
forceVector=cc(:,4);


Load1=[1:numberElements]'; % sorting out the loads of all members in one matrix ""named Load1""
Load=loadmembers; % reading the actual applied loads from excel that already been reed in FrameAnalysis2D

for ii=1:length(Load)    % lood to match and fill the matrix
    c=Load(ii,1);
    Load1(c,2)=Load(ii,2);
end
G = E/(2*(1+V));% shear modulus
L0=memberlength(:,1);
epsilon=sqrt(235/steelGrade);
save('Excel_variables.mat')