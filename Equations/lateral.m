function dxdt = lateral(phi,vphi,r,vr,teta,vteta,aphi,k,Mt,Lc,Tbit,e,ch,...
                                                               Ftheta,Fr)
% lateral Calculates the time derivative of the lateral coordinate vector 
%         for a lateral-torsional lumped parameter model of a drill-string.
% 
%   dxdt = lateral(phi,vphi,r,vr,teta,vteta,aphi,k,Mt,Lc,Tbit,e,ch,...
%                                                                Ftheta,Fr)
%   with dxdt as the time derivative of the lateral coordinate vector.
% 
%   Inputs:
%   phi    -> Torsion angle
%   vphi   -> Drill-bit speed
%   r      -> Radial displacement
%   vr     -> Radial velocity
%   theta  -> Precession angle
%   vtheta -> Precession angular velocity
%   aphi   -> Initial acceleration
%   k      -> Stiffness (N/m)
%   Mt     -> Total mass (kg)
%   Lc     -> BHA section length
%   Tbit   -> Torque on bit
%   e      -> Eccentricity
%   ch     -> Fluid drag damping
%   Ftheta -> Friction force
%   Fr     -> Impact force
%
%   Outputs:
%   dxdt -> Time derivative of the lateral coordinate vector
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


Vt = sqrt(vr^2+(r*vteta)^2);

kl = (k-Tbit*pi^3/(2*Lc^2));   
    

    dxdt(1,1) = vr;
    dxdt(2,1) = -ch*Vt*vr/Mt ...                 % Fluid damping
               + (vteta^2-kl/Mt)*r ...           % Gravity and rigidity
               + e*(vphi^2*cos(phi-teta)  ...    % Eccentricity
               + aphi*sin(phi-teta)) ...         % Eccentricity
               - Fr/Mt;                          % Normal radial force
    
    %vteta
    dxdt(3,1) = vteta;
    
    %ateta
    dxdt(4,1) = -(2*vr/r + ch*Vt/Mt)*vteta ...   % Damping and elasticity               
                + e*(vphi^2*sin(phi-teta) ...    % Eccentricity
                - aphi*cos(phi-teta))/r ...      % Eccentricity
                - Ftheta/r/Mt;                   % Friction force    
end