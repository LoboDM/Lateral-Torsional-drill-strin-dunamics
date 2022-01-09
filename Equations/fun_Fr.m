function Fr = fun_Fr(r,vr,cs,ks,H,tol)     
%  fun_Fr  Calculate the impact force for a laterla-torsional lumped 
%          parameter model of a drill-string.
% 
%   Fr = fun_Fr(r,vr,cs,ks,H,tol)     
%   with Fr as the impact force.
% 
%   Inputs:
%   r   -> Radial displacement
%   vr  -> Radial velocity
%   cs  -> Borehole wall damping
%   ks  -> Borehole wall stiffness
%   H   -> Impact boolean
%   tol -> Gap between collar and borehole wall
%
%   Outputs:
%   Fr -> Impact force
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


    Fr = H*(ks*(r - tol) + cs*vr);
end