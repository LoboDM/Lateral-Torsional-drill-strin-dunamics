function dphidt = torcional(t,phi,vphi,kt,ct,k1,c1,Omega,Tbit,Im,invIm,Tlat)

% torcional   Calculates the time derivative of the torsion angle for a   
%             lateral-torsional lumped parameter model of a drill-string.
% 
%   dphidt = torcional(t,phi,vphi,kt,ct,Omega,Tbit,Im,Tlat)
%   with dphidt as the time derivative of the torsion angle.
% 
%   Inputs:
%   t      -> Initial time
%   phi    -> Torsion angle
%   vphi   -> Drill-bit speed
%   kt     -> Torsional stiffness of each ktp element
%   ct     -> Proportional Damping
%   Omega  -> Rotational speed in rad/s
%   Tbit   -> Torque on bit
%   Im     -> Inertia matrix
%   Tlat   -> Torque due to coupling
%
%   Outputs:
%   dphidt -> Time derivative of the torsion angle
% 
%  LAST MODIFIED: 28/11/2021 BY DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS

n = length(kt);

T = zeros(n,1);
T(1,1) = T(1,1) + c1*Omega + k1*Omega*t;
T(end,1) = T(end,1) - Tbit(end) + Tlat;

% State Matrix
A = [zeros(n) eye(n); -invIm*kt -invIm*ct];

f = [zeros(n,1); invIm*T];

dphidt = A*[phi; vphi] + f;
end