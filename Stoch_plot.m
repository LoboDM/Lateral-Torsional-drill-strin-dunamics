clear all

% INPUT
N_sim      = 250;   % Number of simulations
WOB_plot   = 205;   % WOB simulated
rpm_plot   = 85;   % RPM simulated
bha_region = 1;     % region in BHA to assess

% OUTPUT
mw_spec = 128;  % Window for spectrogram
ocerlap = 120;  % Overlap of windows in spectrogram
% IMPORTANT! The values above are scaled to fit simulation time step

addpath(strcat(pwd,'\app_utils'));

% Files (input and output)
local = strcat(pwd, '\Results_stoch', '\WOB =',num2str(round(WOB_plot)*1000),...
    'rpm =',num2str(rpm_plot,'%03.f'));
output_file = strcat(local,'\WOB =',num2str(round(WOB_plot)*1000),...
    'rpm =',num2str(rpm_plot,'%03.f'),"_stochplot");

input_file = strcat(local,'\WOB =',num2str(round(WOB_plot)*1000),...
    'rpm =',num2str(rpm_plot,'%03.f'),"_",num2str(1));
file_values = load(input_file);

%% Spectrogram parameters calculation
dt = (file_values.t(2)-file_values.t(1));
sp_window = round(mw_spec*1e-2/dt);
sp_noverlap = round(ocerlap*1e-2/dt);
sp_nfft = round(mw_spec*1e-2/dt);

x = file_values.r(bha_region,:).*cos(file_values.theta(bha_region,:));
y = file_values.r(bha_region,:).*sin(file_values.theta(bha_region,:));   

[~,~,t1,~] = spectrogram(x+y*sqrt(-1),sp_window,sp_noverlap,sp_nfft,...
        1/dt,'centered','power','yaxis');

p1_list = zeros(length(t1),N_sim);
p2_list = zeros(length(t1),N_sim);
f_list  = zeros(length(t1),N_sim);
f_list1  = zeros(length(t1),N_sim);
f_list2  = zeros(length(t1),N_sim);

ppm = ParforProgressbar(N_sim, 'progressBarUpdatePeriod', 2.5);
parfor ii_sim2 = 1:N_sim
    disp(['Getting result ' num2str(ii_sim2) '/' num2str(N_sim) '...'])
    
    input_file = strcat(local,'\WOB =',num2str(round(WOB_plot)*1000),'rpm =',...
            num2str(rpm_plot,'%03.f'),"_",num2str(ii_sim2));
    file_values = load(input_file);
    
    x = file_values.r(bha_region,:).*cos(file_values.theta(bha_region,:));
    y = file_values.r(bha_region,:).*sin(file_values.theta(bha_region,:));   

    [~,f,t1,p] = spectrogram(x+y*sqrt(-1),sp_window,sp_noverlap,sp_nfft,...
        1/dt,'centered','power','yaxis');
    
    % Dominant frequency
    [fridge,~,~] = tfridge(p,f);
    
    % Power of dominant positive frequency
    f1 = f > 0;
    p1 = p(f1,:);
    [fridge1,~,lr1] = tfridge(p1,f(f1));
    
    % Power of dominant negative frequency
    f2 = f < 0;
    p2 = p(f2,:);
    [fridge2,~,lr2] = tfridge(p2,f(f2));
    
    f_list(:,ii_sim2) = fridge;
    f_list1(:,ii_sim2) = fridge1;
    f_list2(:,ii_sim2) = fridge2;
    p1_list(:,ii_sim2) = 10*log10(p1(lr1));
    p2_list(:,ii_sim2) = 10*log10(p2(lr2));
    vphi_list(:,ii_sim2) = file_values.vphi;
    ppm.increment();
end
delete(ppm);
save(output_file);
%%
local = strcat(pwd, '\Results_stoch', '\WOB =',num2str(round(WOB_plot)*1000),...
    'rpm =',num2str(rpm_plot,'%03.f'));
output_file = strcat(local,'\WOB =',num2str(round(WOB_plot)*1000),...
    'rpm =',num2str(rpm_plot,'%03.f'),"_stochplot");
load(output_file)

figure(1)
plot(t1,p1_list)
title('f>0')
xlabel('Time (s)')
ylabel('Power (dB)')

figure(2)
plot(t1,p2_list)
title('f<0')
xlabel('Time (s)')
ylabel('Power (dB)')

figure(3)
plot(t1,f_list)
title('Dominant frequency')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

figure(31)
plot(t1,quantile(f_list',0.05)), hold on
plot(t1,mean(f_list'))
plot(t1,quantile(f_list',0.95)), hold off
title('Dominant frequency')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%%
figure(32)
plot(file_values.t,vphi_list*60/2/pi)
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$\dot{\phi}$ (rpm)','Interpreter','latex','FontSize',16)

%%
figure(4)
plot(t1,f_list1)
title('Dominant positive frequency')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

figure(5)
plot(t1,f_list2)
title('Dominant negative frequency')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

figure(6)
p_ratio = log10((10.^(p1_list))./(10.^p2_list));
plot(t1,p_ratio)
title('Magnitude ratio')
xlabel('Time (s)')
ylabel('f>0/f<0 - 1')

figure(7)
f_ratio = abs(f_list1)./abs(f_list2);
plot(t1,f_ratio)
title('Frequency ratio')
xlabel('Time (s)')
ylabel('f>0/f<0')

%% Plots a specific realization ii_sim2

% INPUT
ii_sim2 = 11; % Realization to be plotted

% OUTPUT
local = strcat(pwd, '\Results_stoch', '\WOB =',num2str(round(WOB_plot)*1000),...
    'rpm =',num2str(rpm_plot,'%03.f'));
input_file = strcat(local,'\WOB =',num2str(round(WOB_plot)*1000),'rpm =',...
        num2str(rpm_plot,'%03.f'),"_",num2str(ii_sim2));
load(input_file);
if ~exist(input_file, 'dir')
   mkdir(input_file)
end

addpath(strcat(pwd,'\Equations'));
addpath(strcat(pwd,'\Solver'));
addpath(strcat(pwd,'\Plots'));

% Plot 1 - radial disp, torsional and whirl speeds vs time 
figure(100)
subplot(3,1,1)
plot(t,r(bha_region,:))
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$r$ (m)','Interpreter','latex','FontSize',16)

subplot(3,1,2)
plot(t,vphi*60/2/pi)
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$\dot{\phi}$ (rpm)','Interpreter','latex','FontSize',16)

subplot(3,1,3)
plot(t,vtheta(bha_region,:))
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$\dot{\theta}$ (rad/s)','Interpreter','latex','FontSize',16)
saveas(gca,strcat(input_file, '\radial_disp_tors_whirl_speed'),'fig')

% Plot 2 - Phase diagram
t222 = 40;
x = r(bha_region,t>t222).*cos(theta(bha_region,t>t222));
y = r(bha_region,t>t222).*sin(theta(bha_region,t>t222));   
figure(102)
plot(x,y)
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
saveas(gca,strcat(input_file, '\Phase_diagram'),'fig')

% Plot 3 - Lateral natural frequency vs time
i_dof = LATERAL_dofs(bha_region);
figure(103)
plot(t,(WOB(i_dof)*pi^2/(2*Lc(i_dof)))*ones(length(t),1)), hold on
plot(t,fun_Tbit(vphi,WOBf)*pi^3/(2*Lc(i_dof)^2))
plot(t,(pi^4)*E*I_area(i_dof)./(2*Lc(i_dof).^3)*ones(length(t),1))
plot(t,k(3) - fun_Tbit(vphi,WOBf)*pi^3/(2*Lc(i_dof)^2)), hold off
legend('WOB','TOB','No-load','Total')
temp = k(3) - fun_Tbit(vphi,WOBf)*pi^3/(2*Lc(i_dof)^2);
disp('Lateral natural frequency'), disp(sqrt(temp(end)/Mt(i_dof))/2/pi)
saveas(gca,strcat(input_file, '\Lat_nat_freq_time'),'fig')

% Plot 5 - Spectrogram
x = r(bha_region,:).*cos(theta(bha_region,:));
y = r(bha_region,:).*sin(theta(bha_region,:));   
mw_spec = 128;
ocerlap = 120;
figure(105)
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
saveas(gca,strcat(input_file, '\Spectrogram'),'fig')

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
saveas(gca,strcat(input_file, '\Spectrogram_2'),'fig')

figure(52)
plot(t1,fridge1,'LineWidth',2)

hold on
plot(t1,fridge2,'LineWidth',2)
legend('f>0','f<0')
xlabel('Time (s)')
ylabel('Power (dB)')
hold off
saveas(gca,strcat(input_file, '\Spectrogram_3'),'fig')

% Plot 6 - Bit-rock interaction
dphi_range = 0:0.1:30;
Tbit_range = fun_Tbit(dphi_range,WOBf);

figure(106)
plot(dphi_range*180/pi,Tbit_range/1000)
xlabel('Bit Speed [RPM]')
ylabel('TOB [kN.m]')
saveas(gca,strcat(input_file, '\BitRock_interaction'),'fig')

% Plot 7 - ROP
ROP = fun_ROP(vphi,WOBf);

figure (107)
plot(t,ROP)
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$\dot{z}$ (m/s)','Interpreter','latex','FontSize',16)
saveas(gca,strcat(input_file, '\ROP'),'fig')

% Plot 8 - Axial displacement
figure (108)
plot(t,cumsum(ROP*dt)), hold on
plot(t,z,'--'), hold off
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$z$ (m)','Interpreter','latex','FontSize',16)
saveas(gca,strcat(input_file, '\Axial_Disp'),'fig')

% Plot 9 - Well radius
H_s = Hs_extract(H_grid,z_grid,theta_grid,theta,z);
figure (109)
plot(t,(1+H_s)*Dbwall/2)
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$R_{wall}$','Interpreter','latex','FontSize',16)
saveas(gca,strcat(input_file, '\Well_Radius'),'fig')
%%
t222 = 100;
x = r(bha_region,t>t222).*cos(theta(bha_region,t>t222));
y = r(bha_region,t>t222).*sin(theta(bha_region,t>t222));
t_d = t(t>t222);
vphi_d = vphi(t>t222);
z_d = z(t>t222);

theta_d = linspace(0,2*pi,51);

figure('Position', [200 200 900 300])
subplot(1,2,1)
xlabel('Time (s)')
ylabel('Bit RPM')
gg = animatedline('color','b');
axis([t_d(1), t_d(end), 0, 1.2*max(vphi_d*60/pi)])
grid on
subplot(1,2,2)
g = animatedline('color','b','marker','o');
g2 = animatedline('color','r');

x_wall_ref = tol(LATERAL_dofs(bha_region)).*cos(theta_d);
y_wall_ref = tol(LATERAL_dofs(bha_region)).*sin(theta_d);
hold on, plot(x_wall_ref,y_wall_ref,'k--')

axis square;
axis(1e-3*[-50, 50, -50, 50])
xticks(-0.05:0.025:0.05)
yticks(-0.05:0.025:0.05)
xlabel('X (m)')
ylabel('Y (m)')
grid on
nnn = 50;
lim = length(x)/nnn;

for ii = 1:5:length(x)-nnn
    clearpoints(g)
    addpoints(g,x(ii:ii+nnn),y(ii:ii+nnn))
    addpoints(gg,t_d(ii:ii+nnn),vphi_d(ii:ii+nnn)*60/pi)
    
    H_sd = Hs_extract(H_grid,z_grid,theta_grid,theta_d,...
                                          z_d(ii+nnn)*ones(size(theta_d)));

    clearpoints(g2)
    x_wall = (tol(LATERAL_dofs(bha_region)) + H_sd*Dbwall/2).*cos(theta_d);
    y_wall = (tol(LATERAL_dofs(bha_region)) + H_sd*Dbwall/2).*sin(theta_d);
    addpoints(g2,x_wall,y_wall)
    
    
    drawnow update
    pause(0.03)
end
