function [K,k1] = K_Gen(NgdlTotDP,NgdlTotBHA,kDP,kBHA,i_print)

% K_Gen Generates the stiffness matrix and the Stiffness of the first pipe
%       segment for a torsional lumped parameter model of a drill-string
%       considering the BHA as a rigid or flexible body.
%
%   [K,k1] = K_Gen(NgdlTotDP,NgdlTotBHA,kDP,kBHA)
%   with K as the stiffness matrix and k1 as the Stiffness of the first 
%   pipe segment generated for a model with NgdlTotDP degrees of freedom at
%   the drill-pipes, drill-pipes of stiffness kDP, BHA with stiffness kBHA.
%
%   Inputs:
%   NgdlTotDP   -> Total number of degrees of freedom at the drill pipes
%   NgdlTotBHA  -> Total number of degrees of freedom at the drill collars
%                  If NgdlTotBHA == 0, BHA is a rigid body
%   kDP         -> Drill pipes stiffness
%   kBHA        -> Drill collars stiffness
%   i_print     -> If True, print usefull information
%
%   Outputs:
%   K           -> Stiffness matrix
%   k1          -> Stiffness of the first pipe segment
%
%   Ex:
%   NgdlTotDP = 5, NgdlTotBHA = 0, kDP = [K1 K2 K3 K4 K5], kBHA = [K6 K7]
%   [K,k1] = K_Gen(NgdlTotDP,NgdlTotBHA,NgdlTot,kDP,kBHA) 
%   K = [ K1+K2   -K2     0        0       0;
%          -K2   K2+K3   -K3       0       0;
%           0     -K3   K3+K4     -K4      0;
%           0      0     -K4     K4+K5   -K5;
%           0      0      0       -K5     K5]    
%
%   k1 = K1
%
%   IMPORTANT!! To consider more dofs at the BHA, it is recommended to
%   increase the dofs at the DP as well.
%
%   LAST MODIFIED: 27/03/2020 BY JORDAN BARBOZA AND DANIEL LOBO
%   CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


% vector with the stiffness of each segment of the drill string
kvector = vertcat(kDP,kBHA);

if NgdlTotBHA > 0       %(Flexible BHA)
    Ngdl = NgdlTotDP + NgdlTotBHA;
else % NgdlTotBHA = 0   %(Rigid BHA)
    Ngdl = NgdlTotDP;
end

K = zeros(Ngdl);
for i = 1:Ngdl-1
    K(i,i) = kvector(i,1) + kvector(i+1,1);
end
K(Ngdl,Ngdl) = kvector(Ngdl,1);

for i = 2:Ngdl
    K(i-1,i) = -kvector(i,1);
    K(i,i-1) = -kvector(i,1);
end

k1 = kvector(1);

if i_print
    disp('Torsional stiffness matrix'), disp(K)
end
end