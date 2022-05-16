function R = R_theta(dtheta,l_theta,R_well)
% R_theta  Calculate the correlation function between two points in a
%          circunference
% 
%   R = R_theta(dtheta,l_theta,R_well)
%   with R as the correlation
% 
%   Inputs:
%   dtheta   -> Distance between points in rad
%   l_theta  -> Lenght of correlation in rad
%   R_well   -> Well radius (Usually assumed as bit radius)
%
%   Outputs:
%   R  -> Correlation between two points in a circunference
% 
%  LAST MODIFIED: 12/05/2022 BY DANIEL LOBO
%  CREATED BY DANIEL LOBO

diff_theta = 2*R_well^2*(1 - cos(dtheta));

R = exp(-diff_theta./l_theta^2);

end