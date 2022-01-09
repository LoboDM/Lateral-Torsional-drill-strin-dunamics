function Regimen_Maps(local,vecWOB,vecrpm,bha_region)
% Regimen_Maps Program used to generate the stability maps of a lateral-
%              torsional lumped parameter model of a drill-string
%              considering axial force in the stiffness matrix.
%              
%  Inputs:
%  local      -> folder address with simulation data.
%  vecWOB     -> WOB vector.                        
%  vecrpm     -> RPM vector                        
%  bha_region -> BHA region to be analyzed.    
%
%  LAST MODIFIED: 07/06/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


close all

folder = local;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PRELIMINARY CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Opens the calculated data
mmm = length(vecWOB);
nnn = length(vecrpm);

Whirl = zeros(mmm,nnn);
StickSlip = zeros(mmm,nnn);
Impact   = zeros(mmm,nnn);

hh =1 ;
kk =1 ;
for iii = 1:mmm
    for jjj = 1:nnn
            % Open files
    namestring = strcat('WOB = ',num2str(vecWOB(iii), '%06.f'),'rpm = ',...
        num2str(vecrpm(jjj), '%03.f'),'.mat');

    arquivo = strcat(folder,namestring);
    load(arquivo);
    dt = diff(t);
    dt = dt(1);
    
    if sum(isnan(r(bha_region,:))) > 0
        corte  = find(isnan(r(bha_region,:)),1) - 1;
        r      = r(bha_region,1:corte);
        theta  = teta(bha_region,1:corte);
        vtheta = vtheta(bha_region,1:corte);
        vphi   = vphi(bha_region,1:corte);
        t      = t(bha_region,1:corte);
        tf     = t(end);
    end
    
    % Precession Analysis
    Whirl(iii,jjj) = analysewhirl(r,theta,dt,tf);
    
    % Stick-Slip Analysis
    StickSlip(iii,jjj) = analysestickslip(t,vphi,vecrpm(jjj),dt);
    
    % Impact Analysis
    Impact(iii,jjj) = analyseimpact(t,r,dt,WOBf,vecrpm(jjj));    
    end
end

% Stability
[Map_Regimen,map12 ]= Regimen(Whirl,StickSlip,Impact);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vecWOB2 = vecWOB/1000; 

% Position of figures
pos = [100 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot stability
pos = pos*0.9;
RegMap = figure(1);

axesEst= axes('Parent',RegMap);
hold(axesEst,'on');

% plot map
C = Map_Regimen;
surf(vecrpm,vecWOB2,Map_Regimen)
colormap(map12);
axis([min(vecrpm)  max(vecrpm) vecWOB2(end-1) vecWOB2(1) 0 5])
view(0,90)
xlabel('$\Omega$ (rpm)','Interpreter','latex','FontSize',18)
ylabel('$W_{ob}$ (kN)','Interpreter','latex','FontSize',18)
set(gca,'FontSize',16)

set(gcf, 'Position', [pos 500 400])

saveas(figure(1),'map','png');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function whirl = analysewhirl(r,teta,dt,tf)

%  analysewhirl Determines the precession orientation.
% 
%   whirl = analysewhirl(r,teta,dt,tf)
%   with whirl as the precession orientation flag.
%      
%   Inputs:
%   r     -> Radial displacement
%   teta  -> Precession angle
%   dt    -> Maximum step/Saving step
%   tf    -> Final time
%
%   Outputs:
%   whirl -> Precession orientation flag
%
%  LAST MODIFIED: 14/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS

% Takes the last 20 points (24s +-)
ti = round(0.5*tf/dt);
tf = round(tf/dt);
r = r(ti:tf);
teta = teta(ti:tf);

% Extracts frequency
Fs = 1/dt;             % Sample frequency
N = length(r)-1; 
dF = Fs/N ;   
f = (0):dF:(Fs);       % Axis of a normal fft

% Cartesian coordinates
x = r.*cos(teta);
y = r.*sin(teta);
z = x+1i*y;

% Find the precession fft
freqr = fft(z);

% Write according to precession
FREQx = (f*60-Fs*30);
AMPLIy = abs(fftshift(freqr));

% Finds the predominant frequency
[i, j] = max(AMPLIy);
whirl_aux = FREQx(j);
% plot(FREQx,AMPLIy)
% Says the type of precession
% if Precessao < 0
%     whirl = -1;
% elseif Precessao > 0 
%     whirl = 1;
if whirl_aux < 0
    whirl = 1; % "Backward"
elseif whirl_aux >= 0 
    whirl = 0; % "Forward"
end

end

function impact = analyseimpact(tspan,r,dt,WOB,rpm)

%  analyseimpact Determines the Impact/contact occurrence.
% 
%   impact = analyseimpact(tspan,r,dt,WOB,rpm)
%   with impact as the impact/contact occurrence flag.
% 
%   Inputs:
%   tspan   -> Time information
%   r       -> Radial displacement
%   dt      -> Maximum step/Saving step
%   WOB     -> Weight on Bit in [N]
%   rpm     -> Rotational speed of rotary table
%
%   Outputs:
%   impact  -> Impact/contact occurrence flag
%
%  LAST MODIFIED: 14/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS

% Get the "permanent regime"
impact = 3;
tf = max(tspan);
ti = round(0.50*tf/dt);
tf = round(tf/dt);

tspan = tspan(ti:tf); 
r = r(ti:tf);

rc =  0.0225; % Same for all analyzes
rmed  = mean(r);
rmax  = max(r);
rmin  = min(r);
ii    = (rmax-rmin)/rc;

if rmed >= rc
    impact = 0; % "Permanent contact"
elseif rmed <= rc && rmax >= rc
    impact = 2; % "Impact"
elseif rmed <= rc && rmax < rc
    impact = 1; % "No contact"
end

if impact == 3
    impact = 1;
%     error = ['WOB error = ', num2str(WOB), ' and rpm = ',num2str(rpm)];  
%       mmm = msgbox(error);
end

end



function stickslip = analysestickslip(tspan,vphi,rpm,dt)

%  analysestickslip is the function that determines the Stick-slip 
%                   occurrence.
% 
%   stickslip = analysestickslip(tspan,vphi,rpm,dt)
%   with stickslip as the stick-slip occurrence flag.
% 
%   Inputs:
%   tspan     -> Time information
%   vphi      -> Drill-bit speed
%   rpm       -> Rotational speed of rotary table
%   dt        -> Maximum step/Saving step
%
%   Outputs:
%   stickslip -> Stick-slip occurrence flag
%
%  LAST MODIFIED: 14/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS

% Get the "permanent regime"
stickslip = 3;
tf = max(tspan);
ti = round(0.50*tf/dt);
tf = round(tf/dt);
rot = rpm;
tspan = tspan(ti:tf); 
vphi = vphi(ti:tf)*60/2/pi;
% plot(tspan,vphi)
vphimax  = max(vphi);
vphimin  = min(vphi);
ss    = (vphimax-vphimin)/(2*rot);

if ss >= 0.80
    stickslip = 1; % "Stick-slip"
elseif ss < 1
    stickslip = 0; % "No stick-slip"
end

if stickslip == 3
    stickslip = 1;
%     mmm = msgbox('integration error!');
end   

end



function bluckling = bucklinganalysis(WOB,vphi,k,Lc)

%  buckling Determines the buckling occurrence.
% 
%   bluckling = bucklinganalysis(WOB,vphi,k,Lc)
%   with bluckling as the buckling occurrence flag.
% 
%   Inputs:
%   WOB       -> Weight on Bit in [N]
%   vphi      -> Drill-bit speed
%   k         -> Stiffness in [N/m]
%   Lc        -> BHA section length
%
%   Outputs:
%   bluckling -> Buckling occurrence flag
%
%  LAST MODIFIED: 14/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS

WOBref = 244.2e+3;
muref = 1.4;
b0 = 4.17e+3;
b1 = 1.91;
b2 = 8.5;
b3 = 5.47;
Xt  = b0*muref/WOBref;

vecTbit =  WOB*Xt*(tanh(b1.*vphi) + b2.*vphi./(1+b3.*vphi.^2));
kl = (k-vecTbit*pi^3/(2*Lc^2));

kcrt = min(kl);
if kcrt < 0
    bluckling = 0; % "Buckling"
else
    bluckling = 1; % "No buckling"
end
end
