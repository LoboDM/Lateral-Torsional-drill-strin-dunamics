function R = R_z(dz,l_z)
% R_theta  Calculate the correlation function between two points in z-dir
% 
%   R = R_z(dz,l_z)
%   with R as the correlation
% 
%   Inputs:
%   dz   -> Distance between points in m
%   l_z  -> Lenght of correlation in m
%
%   Outputs:
%   R  -> Correlation between two points in z-direction
% 
%  LAST MODIFIED: 12/05/2022 BY DANIEL LOBO
%  CREATED BY DANIEL LOBO

diff_z = dz.^2;

R = exp(-diff_z/l_z^2);

end