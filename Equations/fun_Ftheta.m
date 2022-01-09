function Ftheta = fun_Ftheta(vphi,r,vtheta,u,Fn,R)      
%  fun_Ftheta Calculate the friction force for a torsional lumped parameter
%             model of a drill-string.
% 
%   Ftheta = fun_Ftheta(vphi,r,vtheta,u,Fn,R)      
%   with Ftheta as the friction force.
% 
%   Inputs:
%   vphi   -> Drill-bit speed
%   r      -> Radial displacement
%   vtheta -> Precession angular velocity
%   u      -> Wall friction coefficient
%   Fn     -> Impact normal force
%   R      -> Collar Radius
%
%   Outputs:
%   Ftheta -> Friction force
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


% Friction force
Ftheta = tanh((vphi*R+r*vtheta)*1e8)*u*Fn;      
end