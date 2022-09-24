clc
clear all
close all

% Stochastic Parameters
R_well = 0.216/2;
Rco = 0.171/2;
l_theta = pi/50;    % Damping constant proportional to inertia
l_z     = 0.25;     % Damping constant proportional to inertia
z_lim   = 2.5;   % Damping constant proportional to inertia
N_theta = 4;   % Damping constant proportional to inertia 
N_z     = 3;   % Damping constant proportional to inertia
mu_H     = 0.05;    % Damping constant proportional to inertia
sigma_H  = 0.03; % Damping constant proportional to inertia

%% Analysis of correlation functions
dtheta = linspace(-2*pi,2*pi,101);
dz = linspace(-z_lim,z_lim,101);

figure(1)
plot(dtheta,R_theta(dtheta,l_theta,R_well),'linewidth',1.5)
xlim([-2*pi,2*pi])
ylabel('R(\Delta \theta)')
xlabel('\Delta \theta')

figure(2)
plot(dz, R_z(dz,l_z),'linewidth',1.5)
ylabel('R(\Delta z)')
xlabel('\Delta z')

%% Generate random field
[z_grid, theta_grid, H_grid] = ...
         Random_Field_Gen(N_theta, N_z, mu_H,sigma_H, l_theta,l_z, ...
                          z_lim, R_well);

figure(3)
    s = pcolor(theta_grid,z_grid,H_grid);
    s.FaceColor = 'interp';
    s.EdgeAlpha = 0.1;
    xlabel('\theta')
    ylabel('z')
%     caxis([0 3])
    colorbar
    colormap('jet')


theta_map = theta_grid(ceil(N_z/2),:);
H_map     = H_grid(ceil(N_z/2),:);

x = ((1+H_map)*R_well).*cos(theta_map);
y = ((1+H_map)*R_well).*sin(theta_map);

x_ref = (R_well).*cos(theta_map);
y_ref = (R_well).*sin(theta_map);

figure(4)
plot([x_ref x_ref(1)],[y_ref y_ref(1)]), hold on
plot([x x(1)],[y y(1)]), hold off
xlim([-0.151 0.151])
ylim([-0.151 0.151])
grid on

tol_x = ((1+H_map)*R_well-Rco).*cos(theta_map);
tol_y = ((1+H_map)*R_well-Rco).*sin(theta_map);


tol_x_ref = (R_well-Rco).*cos(theta_map);
tol_y_ref = (R_well-Rco).*sin(theta_map);

figure(5)
plot([tol_x_ref tol_x_ref(1)],[tol_y_ref tol_y_ref(1)]), hold on
plot([tol_x tol_x(1)],[tol_y tol_y(1)]), hold off
xlim([-0.05 0.05])
ylim([-0.05 0.05])
grid on

x = ((1+H_grid)*R_well).*cos(theta_grid);
y = ((1+H_grid)*R_well).*sin(theta_grid);

x_ref = (R_well).*cos(theta_grid);
y_ref = (R_well).*sin(theta_grid);



%%
figure(6)
surf([x_ref x_ref(:,1)]',[y_ref y_ref(:,1)]',[z_grid z_grid(:,1)]')

figure(7)
surf([x x(:,1)]',[y y(:,1)]',[z_grid z_grid(:,1)]')
