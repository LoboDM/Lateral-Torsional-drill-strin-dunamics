%  PRINCIPAL.m  Program used for anayisis and simulation of a
%               lateral-torsional lumped parameter model of a drill-string
%               Conducts both single simulation and map simulations.
%
%  Input variables:
%  Analysis   -> Analysis type (single or map)                  Type: Float
%  __________________ ex: Analysis = 1;
%  Lp        -> Drill pipes length.                             Type: Float
%  __________________ ex: Lp = 5000;
%  Dpo       -> Drill pipes outer diameter.                     Type: Float
%  __________________ ex: Dpo = 0.160;
%  Dpi       -> Drill pipes inner diameter.                     Type: Float
%  __________________ ex: Dpi = 0.119;
%  E         -> Modulus of elasticity.                          Type: Float
%  __________________ ex: E = 200*10^9;
%  G         -> Shear coefficient.                              Type: Float
%  __________________ ex: G = 90*10^9;
%  rho       -> Specific mass.                                  Type: Float
%  __________________ ex: rho = 8000;
%  rho_f     -> Fluid specific mass.                            Type: Float
%  __________________ ex: rho_f = 2000;
%  Ca        -> Added mass coefficient.                         Type: Float
%  __________________ ex: Ca = 1.5;
%  Cd        -> Drag coefficient.                               Type: Float
%  __________________ ex: Cd = 1.2;
%  alpha     -> Drill-string slope.                             Type: Float
%  __________________ ex: alpha = 0.0002;
%  g         -> Acceleration of gravity.                        Type: Float
%  __________________ ex: g = 9.8;
%  Dbit      -> Bit diameter.                                   Type: Float
%  __________________ ex: Dbit = 0.2;
%  ks        -> Borehole wall stiffness.                        Type: Float
%  __________________ ex: ks = 10e+8;
%  cs        -> Borehole wall damping.                          Type: Float
%  __________________ ex: cs = 0;
%  u         -> Wall friction coefficient.                      Type: Float
%  __________________ ex: u = 0.35;
%  Lc        -> BHA section length.                             Type: List
%  __________________ ex: Lc = [171.30 267.10 8.55 19.25];
%  Dco       -> BHA collar outer diameter.                      Type: List
%  __________________ ex: Dco = [0.140 0.171 0.171 0.171];
%  Dci       -> BHA collar inner diameter.                      Type: List
%  __________________ ex: Dci = [0.076 0.071 0.071 0.071];
%  a_c       -> Damping constant proportional to inertia.       Type: Float
%  __________________ ex: a_c = 0.36;
%  b_c       -> Damping constant proportional to stiffness.     Type: Float
%  __________________ ex: b_c = 0.41;
%  ti        -> Initial time.                                   Type: Float
%  __________________ ex: ti = 0;
%  tf        -> Final time.                                     Type: Float
%  __________________ ex: tf = 50;
%  dt        -> Maximum step/Saving step.                       Type: Float
%  __________________ ex: dt = 1e-2;
%  tolerance -> Solver tolerance.                               Type: Float
%  __________________ ex: tolerance = 0.5e-4;
%  WOBrange  -> Weight on Bit in [kN] range.                    Type: List
%  __________________ ex: WOBrange = 220:-10:0;
%  rpmrange  -> Rotational speed of rotary table range.         Type: List
%  __________________ ex: rpmrange = 40:5:140;
%  WOB       -> Weight on Bit in [kN] for a single simulation.  Type: Float
%  __________________ ex: WOB = 220;
%  rpm       -> Rotational speed at top for single simulation.  Type: Float
%  __________________ ex: rpm = 40;
%  bha_region-> BHA region analysed.                            Type: Int
%  __________________ ex: bha_region = 1;
%  LATERAL_dofs-> Indices from Lc considered in lat. dyn.       Type: List
%  __________________ ex: LATERAL_dofs = [3 4];
%  N_tor      -> Number of torsional DOfs.                      Type: Int
%  __________________ ex: N_tor = 1;
%
%  LAST MODIFIED: 28/11/2021 BY DANIEL LOBO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS AND DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Analysis = 3;   % Analysis = 1 -> single simulation
                % Analysis = 2 -> regimen map simulation

plot_modes  = [1 2];% index of the normal modes you want to plot
                    % IMPORTANT!! The lowest normal mode index corresponds 
                    % to 1 and the highest cannot exceed the number of
                    % degrees of freedom of the system.
                    
%%%%%%%%%%%%%%%%%%%% DRLL-STRING GEOMETRY AND PROPERTIES %%%%%%%%%%%%%%%%%%
E      = 220*10^9;  % Modulus of elasticity
G      = 85.3e+9;   % Shear coefficient
rho    = 7800;      % Specific mass

rho_f  = 1500;      % Fluid specific mass
Ca     = 1.7;       % Added mass coefficient
Cd     = 1.0;       % Drag coefficient

alpha  = 1e-7;      % Drill-string slope
g      = 9.81;      % Acceleration of gravity

Dbit   = 0.216;     % Bit diameter

% Borehole wall properties
ks  =  1e+9;       % Borehole wall stiffness
cs  =  0;           % Borehole wall damping
u   =  0.35;        % Wall friction coefficient

Lp     = 4733.6;    % Drill pipes length
Dpo    = 0.140;     % Drill pipes outer diameter
Dpi    = 0.119;     % Drill pipes inner diameter

Lc  = [171.30 267.10 8.55 19.25];   % BHA section length
Dco = [0.140 0.171 0.171 0.171];    % BHA collar outer diameter
Dci = [0.076 0.071 0.071 0.071];    % BHA collar inner diameter

% Setting the dofs   
N_tor        = 1;      % Number of DOFs for torsional dynamics
LATERAL_dofs = [3];   % DOfs with lateral dynamics

% Rayleigh damping coefficients
a_c = 0.35;     % Damping constant proportional to inertia
b_c = 0.06;     % Damping constant proportional to stiffness

% Stochastic Parameters
l_theta = pi/60; % Damping constant proportional to inertia
l_z     = 0.25;  % Damping constant proportional to inertia
z_lim   = 2.5;   % Damping constant proportional to inertia
N_theta = 101;   % Damping constant proportional to inertia 
N_z     = 101;   % Damping constant proportional to inertia
mu_H     = 0;    % Damping constant proportional to inertia
sigma_H  = 0.15; % Damping constant proportional to inertia

%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial time
ti = 0;

% Final time
tf = 200;

% Maximum step/Saving step
dt = 1e-3;

% Solver tolerance
tolerance = 1e-4;%0.5e-3;

% WOB and rpm range
WOBrange   = 25:25:225;       % in kN 10
rpmrange   = 40:10:140;        % in rpm 5
WOB        = 135;
rpm        = 105;

% Map properties
bha_region = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVING THE PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save directory
local = [pwd '\Results'];
if ~exist(local, 'dir')
       mkdir(local)
end

lLc = length(Lc);
Dbwall = Dbit;      % Borehole diameter is the same as bit

addpath(strcat(pwd,'\Equations'));
addpath(strcat(pwd,'\Solver'));
addpath(strcat(pwd,'\Plots'));

% Treats WOB info to N
WOBrange = WOBrange*1000;
WOB      = WOB*1000;

% Stochastic field calculation
if Analysis == 3
    [z_grid, theta_grid, H_grid] = ...
         Random_Field_Gen(N_theta, N_z, mu_H,sigma_H, l_theta,l_z, ...
                          z_lim, Dbit/2);
else
    z_grid = 0;
    theta_grid = 0;
    H_grid = 1;
end

if ismember(Analysis, [1 3])
    tic
    fun_Model(ti,tf,dt,tolerance,rho,rho_f,Dco,Dci,Ca,Dpo,Dpi,...
              Dbwall,Lp,E,G,Cd,g,u,ks,cs,alpha,a_c,b_c,Lc,rpm,...
              WOB,local,LATERAL_dofs,N_tor,true,z_grid,theta_grid,H_grid);
    file = strcat(local,'\WOB =',num2str(round(WOB)),'rpm =',...
        num2str(rpm,'%03.f'));
    load(file)
    if ~exist(file, 'dir')
       mkdir(file)
    end
    toc
    % Plot 1 - radial disp, torsional and whirl speeds vs time 
    figure(1)
    subplot(3,1,1)
    plot(t,r(bha_region,:))
    xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
    ylabel('$r$ (m)','Interpreter','latex','FontSize',16)

    subplot(3,1,2)
    plot(t,vphi)
    xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
    ylabel('$\dot{\phi}$ (rad/s)','Interpreter','latex','FontSize',16)
    
    subplot(3,1,3)
    plot(t,vtheta(bha_region,:))
    xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
    ylabel('$\dot{\theta}$ (rad/s)','Interpreter','latex','FontSize',16)
    saveas(gca,[file '\radial_disp_tors_whirl_speed'],'fig')

    % Plot 2 - Phase diagram
    t222 = 40;
    x = r(bha_region,t>t222).*cos(theta(bha_region,t>t222));
    y = r(bha_region,t>t222).*sin(theta(bha_region,t>t222));   
    figure(2)
    plot(x,y)
    xlabel('X','Interpreter','latex','FontSize',16)
    ylabel('Y','Interpreter','latex','FontSize',16)
    saveas(gca,[file '\Phase_diagram'],'fig')
    
    % Plot 3 - Lateral natural frequency vs time
    i_dof = LATERAL_dofs(bha_region);
    figure(3)
    plot(t,(WOB(i_dof)*pi^2/(2*Lc(i_dof)))*ones(length(t),1)), hold on
    plot(t,fun_Tbit(vphi,WOBf)*pi^3/(2*Lc(i_dof)^2))
    plot(t,(pi^4)*E*I_area(i_dof)./(2*Lc(i_dof).^3)*ones(length(t),1))
    plot(t,k(3) - fun_Tbit(vphi,WOBf)*pi^3/(2*Lc(i_dof)^2)), hold off
    legend('WOB','TOB','No-load','Total')
    temp = k(3) - fun_Tbit(vphi,WOBf)*pi^3/(2*Lc(i_dof)^2);
    disp('Lateral natural frequency'), disp(sqrt(temp(end)/Mt(i_dof))/2/pi)
    saveas(gca,[file '\Lat_nat_freq_time'],'fig')
    
    % Plot 4 - Natural shapes
    [Eigen_vec,Eigen_val] = eig(kt,Im);
    for ii = 1:length(Eigen_vec)
        figure(100+ii)
        plot([0; Eigen_vec(:,ii)]./max(Eigen_vec(:,ii)),-cumsum([0; LDPv]))
        title(['Torsional Mode ' num2str(ii) ', Freq = ', ...
                              num2str(Eigen_val(ii,ii)/2/pi,'%.2f') ' Hz'])
        ylabel('Distance from top (m)')
        ylabel('Normalized Magnitude')
        saveas(gca,[file '\Mode_' num2str(ii)],'fig')
    end
        
    
    % Plot 5 - Spectrogram
    x = r(bha_region,:).*cos(theta(bha_region,:));
    y = r(bha_region,:).*sin(theta(bha_region,:));   
    mw_spec = 128;
    ocerlap = 120;
    figure(5)
    [~,f,t1,p] = spectrogram(x+y*sqrt(-1),round(mw_spec*1e-2/(t(2)-t(1))),...
        round(ocerlap*1e-2/(t(2)-t(1))),round(mw_spec*1e-2/(t(2)-t(1))),...
        1/(t(2)-t(1)),'centered','power','yaxis');
    [fridge,~,lr] = tfridge(p,f);
    h = pcolor(t1,f,10*log10(p));
    set(h,'EdgeColor','none')
    hCB = colorbar;
    hCB.Title.String = 'Power (dB)';
    ylim([-30 30])
    hold on
    plot3(t1,fridge,p(lr),'LineWidth',2)
    hold off
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    saveas(gca,[file '\Spectrogram'],'fig')
    

    f1 = f > 0;
    p1 = p(f1,:);
    [fridge1,~,lr1] = tfridge(p1,f(f1));
    
    f2 = f < 0;
    p2 = p(f2,:);
    [fridge2,~,lr2] = tfridge(p2,f(f2));
    
    figure(51)
    plot(t1,10*log10(p1(lr1)),'LineWidth',2)
    
    hold on
    plot(t1,10*log10(p2(lr2)),'LineWidth',2)
    legend('f>0','f<0')
    xlabel('Time (s)')
    ylabel('Power (dB)')
    hold off
    

    saveas(gca,[file '\Spectrogram_2'],'fig')
    
    % Plot 6 - Bit-rock interaction
    dphi_range = 0:0.1:30;
    Tbit_range = fun_Tbit(dphi_range,WOBf);
    
    figure(6)
    plot(dphi_range*180/pi,Tbit_range/1000)
    xlabel('Bit Speed [RPM]')
    ylabel('TOB [kN.m]')
    saveas(gca,[file '\BitRock_interaction'],'fig')
    
    % Plot 7 - ROP
    ROP = fun_ROP(vphi,WOBf);
    
    figure (7)
    plot(t,ROP)
    xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
    ylabel('$\dot{z}$ (m/s)','Interpreter','latex','FontSize',16)
    saveas(gca,[file '\ROP'],'fig')
    
    
    % Plot 8 - Axial displacement
    figure (8)
    plot(t,cumsum(ROP*dt)), hold on
    plot(t,z,'--'), hold off
    xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
    ylabel('$z$ (m)','Interpreter','latex','FontSize',16)
    saveas(gca,[file '\Axial_Disp'],'fig')
          
elseif Analysis == 2
    for i = 1:length(rpmrange)
        parfor j = 1:length(WOBrange)
            disp([num2str((i-1)*length(WOBrange)+j) ' / ' num2str(length(WOBrange)*length(rpmrange))])
            tic
            fun_Model(ti,tf,dt,tolerance,rho,rho_f,Dco,Dci,Ca,Dpo,Dpi,...
                Dbwall,Lp,E,G,Cd,g,u,ks,cs,alpha,a_c,b_c,Lc,rpmrange(i),...
                WOBrange(j),local,LATERAL_dofs,N_tor,false,false);
            toc
        end
    end
    Regimen_Maps(local,WOBrange,rpmrange,bha_region)

end



% figure
% g = animatedline('color','b','marker','o');
% axis(1e-3*[-25, 25, -25, 25])
% nnn = 200;
% lim = length(x)/nnn;
% for ii = 48/dt:5:length(x)-nnn
%     clearpoints(g)
%     addpoints(g,x(ii:ii+nnn),y(ii:ii+nnn))
%     drawnow update
%     pause(0.05)
% end