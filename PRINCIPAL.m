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
Analysis = 1;   % Analysis = 1 -> single simulation
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
ks  =  10e+8;       % Borehole wall stiffness
cs  =  0;           % Borehole wall damping
u   =  0.35;        % Wall friction coefficient

Lp     = 4733.6;    % Drill pipes length
Dpo    = 0.140;     % Drill pipes outer diameter
Dpi    = 0.119;     % Drill pipes inner diameter

Lc  = [171.30 267.10 8.55 19.25];   % BHA section length
Dco = [0.140 0.171 0.171 0.171];    % BHA collar outer diameter
Dci = [0.076 0.071 0.071 0.071];    % BHA collar inner diameter

% Setting the dofs
N_tor        = 2;      % Number of DOFs for torsional dynamics
LATERAL_dofs = [3 4];   % DOfs with lateral dynamics

% Rayleigh damping coefficients
a_c = 0.35;     % Damping constant proportional to inertia
b_c = 0.06;     % Damping constant proportional to stiffness

%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial time
ti = 0;

% Final time
tf = 50;

% Maximum step/Saving step
dt = 1e-2;

% Solver tolerance
tolerance = 1e-4;%0.5e-3;

% WOB and rpm range
WOBrange   = 220:-10:190;       % in kN
rpmrange   = 120:10:140;        % in rpm
WOB        = 200;
rpm        = 90;

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

if Analysis == 1
    fun_Model(ti,tf,dt,tolerance,rho,rho_f,Dco,Dci,Ca,Dpo,Dpi,...
              Dbwall,Lp,E,G,Cd,g,u,ks,cs,alpha,a_c,b_c,Lc,rpm,...
              WOB,local,LATERAL_dofs,N_tor,true);
    file = strcat(local,'\WOB =',num2str(round(WOB)),'rpm =',...
        num2str(rpm,'%03.f'));
    load(file)
    
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
    
    % Plot 2 - Phase diagram
    x = r(bha_region,t>20).*cos(theta(bha_region,t>20));
    y = r(bha_region,t>20).*sin(theta(bha_region,t>20));   
    figure(2)
    plot(x,y)
    xlabel('X','Interpreter','latex','FontSize',16)
    ylabel('Y','Interpreter','latex','FontSize',16)
    
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
    
    % Plot 4 - Natural shapes
    [Eigen_vec,Eigen_val] = eig(kt,Im);
    for ii = 1:length(Eigen_vec)
        figure(100+ii)
        plot([0; Eigen_vec(:,ii)]./max(Eigen_vec(:,ii)),-cumsum([0; LDPv]))
        title(['Torsional Mode ' num2str(ii) ', Freq = ', ...
                              num2str(Eigen_val(ii,ii)/2/pi,'%.2f') ' Hz'])
        ylabel('Distance from top (m)')
        ylabel('Normalized Magnitude')
    end
        
    
    % Plot 5 - Spectrogram
    x = r(bha_region,:).*cos(theta(bha_region,:));
    y = r(bha_region,:).*sin(theta(bha_region,:));   
    figure(5)
    spectrogram(x+y*sqrt(-1),128,120,128,1/(t(2)-t(1)),'centered','power','yaxis')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate nat frequencies, associated eigenvectors and normal modes
    [V,w,MN] = FreqNat_Gen(kt,Im,NgdlTotDP,0,1);

    L = sum([Lp Lc]);
    if NgdlTotDP == 1
        Coord = [0;LDP;L];
        Coord_inv = L*ones(length(Coord),1) - Coord;
    else
        Coord = zeros(NgdlTotDP+2,1);
        for i = 1:NgdlTotDP
            Coord(i+1) = Coord(i) + LDPv(i);
        end
        Coord(end) = L;
        Coord_inv = L*ones(length(Coord),1) - Coord;
    end
    
    for i = 1:length(plot_modes)
        figure(100+i)
        if NgdlTotDP > 1
            a = [ 0 ; MN(:,plot_modes(i)); MN(end,plot_modes(i))];
            plot(a,Coord_inv,'ks--')
            ylabel('Distance from bit (m)','FontSize',20)
            xlabel('Amplitude','FontSize',20)
            xlim([ -2.0 2.0] )
            ylim([ 0 5500 ])
            set(gca,'fontsize',14)
            title([num2str(NgdlTotDP) ' dofs, '...
                num2str(plot_modes(i)) '^o mode: '...
                num2str(w(plot_modes(i))) ' Hz'],'FontSize',16)
            legend([num2str(plot_modes(i)) '^o mode: '...
                num2str(w(plot_modes(i)))...
                ' Hz'],'Location','northwest')
%             if i_save == 1 || i_save == 3
%                 saveas(gca,['temp/Rigid_BHA_' num2str(NgdlDP)...
%                     'DP_Modo' num2str(plot_modes(i))],'epsc')
%                 saveas(gca,['temp/Rigid_BHA_' num2str(NgdlDP)...
%                     'DP_Modo' num2str(plot_modes(i))],'fig')
%                 saveas(gca,['temp/Rigid_BHA_' num2str(NgdlDP)...
%                     'DP_Modo' num2str(plot_modes(i))],'jpg')
%             end
        else
            a = [ 0 ; MN; MN];
            plot(a,Coord_inv,'ks--')
            ylabel('Distance from bit (m)','FontSize',20)
            xlabel('Amplitude','FontSize',20)
            xlim([ -2.0 2.0] )
            ylim([ 0 5500 ])
            set(gca,'fontsize',14)
            title([num2str(NgdlTotDP)...
                ' dofs, ' num2str(1) '^o mode: ' num2str(w) ' Hz'],...
                'FontSize',16)
            legend([num2str(plot_modes(i)) '^o mode: '...
                num2str(w) ' Hz'],'Location','northwest')
%             if i_save == 1 || i_save == 3
%                 saveas(gca,'temp/KEeq_Rigido','epsc')
%                 saveas(gca,'temp/KEeq_Rigido','fig')
%                 saveas(gca,'temp/KEeq_Rigido','jpg')
%             end
        end
    
    end
        
        
    
    

    %%
elseif Analysis == 2
    for i = 1:length(rpmrange)
        for j = 1:length(WOBrange)
            fun_Model(ti,tf,dt,tolerance,rho,rho_f,Dco,Dci,Ca,Dpo,Dpi,...
                Dbwall,Lp,E,G,Cd,g,u,ks,cs,alpha,a_c,b_c,Lc,rpmrange(i),...
                WOBrange(j),local,LATERAL_dofs,N_tor,false);
        end
    end
    Regimen_Maps(local,WOBrange,rpmrange,bha_region)
end