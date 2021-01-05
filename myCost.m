function OptimumDesign = myCost(UCs,UBs,x,params,Type)
% global numberElement      s elementNodes SectionType ll
% Type determines variable optimized and printing:
% Has formula XY, where X: {0,1} (do not, or show costs) and Y:{1,2,3} (weight, energy, both)

w=0;
e_st=0;
e_pa=0;
e_tt=0;
e_ee=0;

e_cs=0;
e_cp=0;
e_ct=0;
e_ce=0;

Cm=0;
Ct=0;
Cf=0;
Cc=0;
Ce=0;

EE=0;
EC=0;
Cost=0;

numberElements = params.numberElements;
elementNodes = params.elementNodes;

ll = params.ll;
p = 7.850; % [tonne/m3]

for i=1:numberElements
    SectionType=elementNodes(:,3);
    sc = SectionType (i,1);
    
    if sc==2
        choice_column = x(i);  

	% Weight part
     mass=UCs(choice_column,1); % read
    A=UCs(choice_column,27); % read  % just added.........
	a=UCs(choice_column,13); % read  % just added.........
	
	% Energy part
	A=UCs(choice_column,27); % read
	a=UCs(choice_column,13); % read
    mass=UCs(choice_column,1); % read   % just added.........

    else
        choice_beam = x(i);
        mass=UBs(choice_beam,1);
    A=UBs(choice_beam,27); % read  % just added.........
	a=UBs(choice_beam,13); % read  % just added.........
%     % Energy part
% 	A=UBs(choice_column,27); % read
% 	a=UBs(choice_column,13); % read
    
	% Energy part
	A=UBs(choice_beam,27); % read
	a=UBs(choice_beam,13); % read
    mass=UBs(choice_beam,1);   % just added...........

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

   
    EC=e_cs+e_cp+e_ct+e_ce;
    
    % cost
    Cm = Cm+900*mass*L; % S275>1550 AND S355>1638
    Ct = Ct+50*mass*L; % in Economic page ...
    Ce = Ce+243*mass*L; % in spon's page 326
    Cf = Cf + 16*a*L; % in spon's page 172-38
    Cc = Cc + 8.9*a*L; % in spon's page 171
    
    Cost=Cm+Ct+Ce+Cf+Cc;
    
end
switch rem(Type,10) % switch last digit
case 1
	OptimumDesign = w;
case 2
	OptimumDesign = EE;
case 3
	OptimumDesign = [w EE];
otherwise
	OptimumDesign = 1;
end

if Type>9
    fprintf('\n\n-----------------------------------------------------------------------\n');
    disp(['Candidate: ' num2str(x)])
    fprintf('\nWeight: %2.2f\nEmbodied energy: %2.3fx10^3\nEmbodied carbon: %2.3fx10^3\nCost: %4.1f',w,EE/1e3,EC/1e3,Cost);
    fprintf('\n--------------------------------------------------------\n');
end
    
end

