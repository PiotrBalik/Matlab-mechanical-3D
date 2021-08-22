function [MRx, FBRPartY, FBRPartZ,LTBRC1,LTBRC2,CRc,SRx,Defl]  = parametersDesign(Ux,i,sc,id,DesignParams)
% Computes material parameters for given beam/column
% Ux - table of values for either UCs or UBs
% i - the row of data
% sc - section type
% id - element number

% parameter intake
elementNodes = DesignParams.elementNodes;
ggg = DesignParams.ggg;   Ned = DesignParams.Ned;
gammaM0 = DesignParams.gammaM0; gammaM1 = DesignParams.gammaM1;
steelGrade = DesignParams.steelGrade;
epsilon = DesignParams.epsilon;
ita = DesignParams.ita; Ved = DesignParams.Ved;
MYed = DesignParams.MYed;   MZed = DesignParams.MYed;
L = DesignParams.L; C1 = DesignParams.C1;
E = DesignParams.E; G = DesignParams.G;
lambdaDashLT0 = DesignParams.lambdaDashLT0;
moment1 = DesignParams.moment1;
moment3 = DesignParams.moment3;
moment5 = DesignParams.moment5;
beta = DesignParams.beta;
% phiColumnYY = DesignParams.phiColumnYY;
% phiColumnZZ = DesignParams.phiColumnZZ;
% phiLT = DesignParams.phiLT;
% psai = DesignParams.psai;
% psi = DesignParams.psi;
LiveL = DesignParams.LiveL;

tw=Ux(i,4);	h=Ux(i,2);
tf=Ux(i,5); b=Ux(i,3);
cwtw=Ux(i,8);   r=Ux(i,6);
cftf=Ux(i,9);	d=Ux(i,7);
Iw=Ux(i,25);	It=Ux(i,26);
Iy=Ux(i,15);	iy=Ux(i,17);
Iz=Ux(i,16);	iz=Ux(i,18);
Wely=Ux(i,19);	Welz=Ux(i,20);
Wply=Ux(i,21);	Wplz=Ux(i,22);
A=Ux(i,27);

%--------------yeildstrength reduction fy--------------------
tff=tf*1000; % to get it back to normal unit mm

if steelGrade==275
    
    fy(tff(:,1)<=16,1)=275;
    fy(tff(:,1)>16,1)=265;
    fy(tff(:,1)>40,1)=255;
    fy(tff(:,1)>63,1)=245;
    
else     % SteelGrade==355
    fy(tff(:,1)<=16,2)=355;
    fy(tff(:,1)>16,2)=345;
    fy(tff(:,1)>40,2)=335;
    fy(tff(:,1)>63,2)=325;
    fy=fy(:,2);
end

fys=fy*1000; % to convert n/mm2 to kN/m2 , because all the units are meter and kN.

%--------------section classification--------------------

alfa=0.5*(1+(Ned/(fys*d*tw)));  % alfa equations
if alfa >1 || alfa <=-1
    alfa=1;
end

psi=((2*Ned)/(A*fys))-1;   % psi equations

%-----------                   -----------------
sectionslassification=3;
if cftf<=(9*epsilon)
    if alfa>0.5 && cwtw<=((396*epsilon)/(13*alfa-1))
        sectionclassification = 1;
    elseif alfa<=0.5 && cwtw<=((36*epsilon)/(alfa))
        sectionclassification = 1;
    end
    
    
elseif cftf<=(10*epsilon) && cftf>(9*epsilon)
    if alfa>0.5 && cwtw<=((456*epsilon)/(13*alfa-1))
        sectionclassification = 2;
    elseif alfa<=0.5 && cwfw<=((41.5*epsilon)/(alfa))
        sectionclassification = 2;
    end
    
elseif cftf<=(14*epsilon) && cftf>(10*epsilon)
    if psi>-1 && cwtw<=((42*epsilon)/(0.67+0.33*psi))
        sectionclassification = 3;
    elseif psi<=-1 && cwtw<=((62*epsilon)*(1-psi)*(sqrt(-psi)))
        sectionclassification = 3;
    end
end

%--------------cross-sectional resistance--------------------
% $$$ compression ressistance

NcRd=(A*fys)/(gammaM0);

CRc=Ned/NcRd; %<=1  CRc  is a breviation for compression resistance for both beam and column

Av=A-2*b*tf+(tw+2*r)*tf;
hw=h-2*tf;
if Av<ita*hw*tw
    Av=ita*hw*tw;
end

VcRd=(Av*(fys/sqrt(3)))/gammaM0;
SRx=Ved/VcRd; %<=1  SR is a breviation for Shear Ressistance for beam and column
%---------------------------------------------

lambda1=93.9*epsilon;

Mcr=C1*(pi^2*E*Iz/L^2)*sqrt((Iw/Iz)+((L^2*G*It)/(pi^2*E*Iz)));


% flexural buckling about y-y

lambdaDashPartYY=(L/iy)*(1/lambda1); %Lcr needs to be found in case of there is a support in mid of Lcrz or Lcry----------+++++++++++++++++++++++++++++++++++++
            
% alpha   % bucklingcurve -- to get the Alpha imperfection factor.....

if h/b>1.2
    if tf<=0.04
        bucklingcurvePartYY='a';
        bucklingcurvePartZZ='b';
    end
    if tf>0.04
        bucklingcurvePartYY='b';
        bucklingcurvePartZZ='c';
    end
    
elseif h/b<=1.2
    if tf<=0.1
        bucklingcurvePartYY='b';
        bucklingcurvePartZZ='c';
    end
    if tf>0.1
        bucklingcurvePartYY='d';  % There is no c, in 2D structures, it is for buckling about z-z
        bucklingcurvePartZZ='d';  % There is no c, in 2D structures, it is for buckling about z-z
    end
    
end

% alpha -- imperfection factor for beam about Y-Y

if bucklingcurvePartYY=='a'
    alphaPartYY=0.21;
elseif bucklingcurvePartYY=='b'
    alphaPartYY=0.34;
else
    alphaPartYY=0.76;
end
if bucklingcurvePartZZ=='b'
    alphaPartZZ=0.34;
end
if bucklingcurvePartZZ=='c'
    alphaPartZZ=0.49;
end
if bucklingcurvePartZZ=='d'
    alphaPartZZ=0.76;
end
phiPartYY=0.5*(1+alphaPartYY*(lambdaDashPartYY-0.2)+lambdaDashPartYY^2);
chiPartYY=1/(phiPartYY+sqrt(phiPartYY^2-lambdaDashPartYY^2));

NbyRd=chiPartYY*A*fys/gammaM1;

FBRPartY=(Ned/NbyRd); % Flexural buckling Resistance about y-y axis

lambdaDashPartZZ=(L/iz)*(1/lambda1); %Lcr needs to be found----------+++++++++++++++++++++++++++++++++++++

% alpha -- imperfection factor for beam about Z-Z
phiPartZZ=0.5*(1+alphaPartZZ*(lambdaDashPartZZ-0.2)+lambdaDashPartZZ^2);
chiPartZZ=1/(phiPartZZ+sqrt(phiPartZZ^2-lambdaDashPartZZ^2));

NbzRd=chiPartZZ*A*fys/gammaM1;

FBRPartZ=(Ned/NbzRd); % Flexural buckling Resistance about y-y axis

% alphaLT   % bucklingcurve -- to get the AlphaLT imperfection factor.....
% alphaLT -- imperfection factor LT

if h/b<=2
    bucklingcurveLT='b';
    alphaLT=0.21;
end

if h/b>1.2 && h/b<=3.1
    bucklingcurveLT='c';
    alphaLT=0.34;
end

if h/b>3.1
    bucklingcurveLT='d';  % There is no c, in 2D structures, it is for buckling about z-z
    alphaLT=0.76;
end


if abs(moment1)>=abs(moment5)
    psai=moment5/moment1; % psi ,, I added a to psi, to distinguish between psi of section classification and this psi
else
    psai=moment1/moment5;
end

Cmy = 0.4; %default values
if sc==2 % if it's column      % I CHANGEEEEEEEEEEEEED THE SECTIONTYPE >>>>>>>>>>>>>>>>>>>>>>>>>>>
% it wassss >> if sectionslassification==2 % if it's column
    
    if psai>=-1 && psai<=1
        Cmy=0.6+(0.4*psai);
    else
        Cmy=0.4;
    end
    
    if Cmy<0.4
        Cmy=0.4;
    end

    HDD=abs(ggg(elementNodes(id,1))-ggg(elementNodes(id,2)));
    
    sigmaHMax=L/300; %300
    
    sigmaH=HDD;  % DEFLECTION FOR FIXED BEAM WITH UDL, W*L^2/384*E*I
    
    HD=sigmaH/sigmaHMax;
    
else % beam
    if abs(moment1)>=abs(moment5)
        as=moment3/moment1; % as is abreviation for AlphaS
    else
        as=moment3/moment5;
    end

    if as>=0 && as<=1
        if psai>=-1 && psai<=1
            Cmy=0.2+(0.8*as);
        end
    end
    if as>=-1 && as<0
        if psai>=0 && psai<=1
            Cmy=0.1-(0.8*as);
        end
        if psai>=-1 && psai<0
            Cmy=0.1*(1-psai)-(0.8*as);
         end
    end
    
	if Cmy<0.4
         Cmy=0.4;   % THIS WAS NOT ACTIVE IN 2D >>>>>>>>>>>>>>>>>>>>>>>>
    end
    
%     if as <-1
%         Cmy = 0.4;        THIS WAS ACTIVE IN 2D >>>>>>>>>>>>>>>>>>>>>>>>>>
%     end
 
    sigmaVMax=L/360; %360
    sigmaV=((LiveL*L^4)/ (384*E*Iy));  % DEFLECTION FOR FIXED BEAM WITH UDL, W*L^2/384*E*I
    
    VD=sigmaV/sigmaVMax;
end
Cmz=Cmy;  % Cmz = Cmy =CmLT. for columns and beams Acoording to EC3 page 80 ................
%------------------------------------------------
% $$$ Bending Moment  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

RO=(((2*Ved)/(VcRd))-1)^2; %NEW
fyr=((1-RO)*fys);  %NEW

if sectionslassification==1 || sectionslassification==2
    Kyy = Cmy*(1 + ((lambdaDashPartYY-0.2)*(Ned/NbyRd)));
    Kzz=Cmz*(1+(((2*lambdaDashPartZZ)-0.6)*(Ned/NbzRd)));

    if Kyy > Cmy*(1+(0.8*(Ned/NbyRd)))
        Kyy=Cmy*(1+(0.8*(Ned/NbyRd)));
    end
    if Kzz > Cmz*(1+(1.4*(Ned/NbzRd)))
        Kzz=Cmz*(1+(1.4*(Ned/NbzRd)));
    end
    Kyz=0.6*Kzz;
    Kzy=0.6*Kyy;
    
    Wy=Wply;
    McRd=(Wy*fys/gammaM0);
    if Ved>=VcRd/2
        McRd=(Wy*fyr/gammaM0);  % McRd is a moment resistance after the reduction of fy, (due to the presence of shear)
    end
     McbzRd=(Wplz*fys/gammaM1);  % calculate the Mcb,z,Rd in case of class 1 & 2 ............
     
else  % when sectionslassification==3
    Kyy = Cmy*(1 + (0.6*lambdaDashPartYY*(Ned/NbyRd)));
    Kzz=Cmz*(1+(0.6*lambdaDashPartZZ*(Ned/NbzRd)));

    if Kyy > (1+(0.6*(Ned/NbyRd)))
        Kyy = (1+(0.6*(Ned/NbyRd)));
    end
    if Kzz > (1+(0.6*(Ned/NbzRd)))
       Kzz = (1+(0.6*(Ned/NbzRd)));
    end
    Kyz=Kzz;
    Kzy=0.8*Kyy;

    Wy=Wely;
    McRd=(Wy*fys/gammaM0);
    if Ved>=VcRd/2    
        McRd=(Wy*fyr/gammaM0);  % McRd is a moment resistance after the reduction of fy, (due to the presence of shear)
    end
    McbzRd=(Welz*fys/gammaM1); % calculate the Mcb,z,Rd in case of class 3..........
end
lambdaDashLT=sqrt(Wy*fys/Mcr);

phiLT=0.5*(1+alphaLT*(lambdaDashLT-lambdaDashLT0)+beta*lambdaDashLT^2);

chiLT=1/(phiLT+sqrt(phiLT^2-beta*lambdaDashLT^2)); % must be <=1

if chiLT>1
    chiLT=1;
end
% MbRd - Lateral torsional buckling **************************************
MbRd=chiLT*Wy*fys/gammaM1;

MRx=MZed/McRd; %<=1  MR is a breviation for Moment resistance

LTBRC1=(Ned/NbyRd)+(Kyy*(MYed/MbRd))+(Kyz*(MZed/McbzRd));
LTBRC2=(Ned/NbzRd)+(Kzy*(MYed/MbRd))+(Kzz*(MZed/McbzRd));

if sc == 2 %column
    Defl = HD;
else
    Defl = VD;
end

end