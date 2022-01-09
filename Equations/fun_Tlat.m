function Tlat = fun_Tlat(phi,r,vr,theta,vtheta,ch,e,R,Fn,Ft)
%  fun_Tlat Calculates lateral coupling torques for a torsional lumped 
%           parameter model of a drill-string.
% 
%   Tlat = fun_Tlat(phi,r,vr,theta,vtheta,ch,e,R,Im,Fn,Ft)
%   with Tlat as the torque due to lateral-torsional coupling
% 
%   Inputs:
%   phi    -> Rotational angle
%   r      -> Radial displacement
%   vr     -> Radial velocity
%   theta  -> Whirl angle
%   vtheta -> Whirl velocity
%   ch     -> Fluid drag damping
%   e      -> Eccentricity
%   R      -> Collar Radius
%   Im     -> Inertia matrix
%   Fn     -> Impact normal force
%   Ft     -> Impact friction force
%
%   Outputs:
%   Tlat   -> Torque due to lateral-torsional coupling
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


% Tangencial velocity
Vt = sqrt(vr^2+(vtheta*r)^2); 
    
Tlat =(-ch*Vt*vr*e*sin(phi-theta)...    % Radial eccentricity component
    +ch*Vt*r*vtheta*e*cos(phi-theta)... % Angular component of eccentricity
        -Ft*(R-e*cos(phi-theta))...
        -Fn*e*sin(phi-theta));       % Radial impact component

end