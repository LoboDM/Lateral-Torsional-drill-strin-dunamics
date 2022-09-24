function [z_grid, theta_grid, H_grid] = ...
         Random_Field_Gen(N_theta, N_z, mu_H,sigma_H, l_theta,l_z, ...
                          z_lim, R_well)
% Random_Field_Gen  Generate a random field
% 
%   H_grid = Random_Field_Gen(N_theta, N_z, mu_H, sigma_H, l_theta,...
%                             l_z, z_lim)
%   with H_grid as the random field at specified nodes
% 
%   Inputs:
%   N_theta   -> Number of nodes in theta-direction
%   N_z       -> Number of nodes in z-direction
%   mu_H      -> Mean of stochastic field H
%   sigma_H   -> Standard deviation of stochastic field H
%   l_theta   -> Lenght of correlation in rad
%   l_z       -> Lenght of correlation in z-direction in m
%   z_lim     -> Upper limit of penetration in theta-direction in m
%   R_well    -> Well radius
%
%   Outputs:
%   H_grid   -> Random field realization at specified nodes
% 
%  LAST MODIFIED: 12/05/2022 BY DANIEL LOBO
%  CREATED BY DANIEL LOBO

mu_h = log(mu_H) - 0.5*log(1 + sigma_H^2/mu_H^2);
Sigma_h = sqrt(log(1 + sigma_H^2/mu_H^2));

theta   = linspace(0,2*pi,N_theta+1);
theta   = theta(1:end-1);
z       = linspace(-0.5,z_lim,N_z);

[theta_grid,z_grid] = meshgrid(theta,z);
theta_map = theta_grid(:);
z_map     = z_grid(:);

% Mount Covariance Matrix
Cov = zeros(N_theta*N_z);
for ii = 1:N_theta*N_z
    for jj = 1:N_theta*N_z
        
        Cov(ii,jj) = Sigma_h^2*R_theta(theta_map(ii)-theta_map(jj),...
            l_theta,R_well)*...
            R_z(z_map(ii)-z_map(jj),l_z);
        
    end
end
% Cov = Sigma_h^2*eye(N_theta*N_z);
Coeff = 1e-12;

Cov = Cov + eye(N_theta*N_z)*Coeff;

U = chol(sparse(Cov));

Phi = mu_h + U'*randn(N_theta*N_z,1);
H = exp(Phi);

H_grid = reshape(H,N_z,N_theta);

end