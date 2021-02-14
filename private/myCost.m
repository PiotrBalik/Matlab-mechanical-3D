function OptimumDesign = myCost(UCs,UBs,x,params,Type)
% Type argument determines variable optimized and printing:
% Has formula XYZ, where
% X: {0,1} (unconstrained or constrained)
% Y: {0,1} (do not, or show costs)
% Z:{1,2,3} (weight, energy, both)
% ie. 13 - show both, 102 - constr. energy

w=0;
e_st=0; e_pa=0;
e_tt=0; e_ee=0;

e_cs=0; e_cp=0;
e_ct=0; e_ce=0;

Cm=0;   Ct=0;
Cf=0;   Cc=0;
Ce=0;

EE=0;   EC=0;
Cost=0;

numberElements = params.numberElements;
elementNodes = params.elementNodes;

ll = params.ll;
p = 7.850; % [tonne/m3]

for i=1:numberElements
    SectionType=elementNodes(:,3);
    sc = SectionType (i,1);
    choice = x(i);  

    
    if sc==2

	A=UCs(choice,27); % read
	a=UCs(choice,13); % read
    mass=UCs(choice,1); % read   % just added.........

    else

	A=UBs(choice,27); % read
	a=UBs(choice,13); % read
    mass=UBs(choice,1);   % just added...........

    end
    L = ll(i);
    
    % weight value
    w = w + mass*L;  % total mass accumulation

    % energy value
    e_st = e_st + (20.1*1000*mass*L); % material
    e_pa = e_pa + 21*a*L; % painting (fire & Corrosion)
    e_tt = e_tt + (1.66*101.85*mass*L); % 0.97*100km*1.05factor   Transportation
    e_ee = e_ee + (958.672*mass*L); % (0.17*163.6*0.54/0.84)*38.3*0.0014*1000
    EE= e_st+e_pa+e_tt+e_ee;

    % carbon
    e_cs = e_cs + 1.6275*1000*mass*L; %1.55*1.05*w;
    e_cp = e_cp + 0.87*a*L; %Carbon for painting
    e_ct = e_ct + 1.66*7.5075*mass*L; % 0.0715*100*1.05 Transportation carbon
    e_ce = e_ce + 72.9397*mass*L; % (0.17*163.6*0.54/0.84)*2.914*0.0014*1000 Erection

   
    EC=e_cs+e_cp+e_ct+e_ce; %UNUSED? for future?
    
    % cost
    Cm = Cm+900*mass*L; % S275>1550 AND S355>1638
    Ct = Ct+50*mass*L; % in Economic page ...
    Ce = Ce+243*mass*L; % in spon's page 326
    Cf = Cf + 16*a*L; % in spon's page 172-38
    Cc = Cc + 8.9*a*L; % in spon's page 171
    
    Cost=Cm+Ct+Ce+Cf+Cc; %UNUSED? for future?
end% for i

flag=0;
if Type>99
[c,~] = constraints(x,params);
flag = any(c >= 1e-4 ); % if there is violation of constraints, with tolerance
end
switch rem(Type,10) % switch last digit
    case 1
        OptimumDesign = w*(1+flag*1e6);
    case 2
        OptimumDesign = EE*(1+flag*1e10);
    case 3
        OptimumDesign = [w EE].*(1+flag*[1e6 1e10]);
    otherwise
        OptimumDesign = 1;
end

if rem(Type,100)>9
    sep=['\n' repmat('--',[1 35]) '\n'];
    fprintf(['\n' sep]);
    disp(['Candidate: ' num2str(x)])
    fprintf('\nWeight: %2.2f\nEmbodied energy: %2.3fx10^3\nEmbodied carbon: %2.3fx10^3\nCost: %4.1f',w,EE/1e3,EC/1e3,Cost);
    fprintf(sep);
end
    
end