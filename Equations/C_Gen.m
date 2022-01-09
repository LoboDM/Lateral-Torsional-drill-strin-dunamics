function [C,c1] = C_Gen(I,K,alfa,beta,k1,i_print)

% C_Gen Generates the damping matrix and the damping of the first pipe
%       segment for a torsional lumped parameter model of a drill-string
%       considering the BHA as a rigid or flexible body.
%
%   [C,c1] = C_Cen(I,K,Xi,NgdlTot,alfa,beta,k1)
%   with C as the damping matrix and c1 as the damping of the first pipe
%   segment generated for a model with inertia matrix I, stiffness matrix
%   K, damping coefficient Xi for the first moment (only used in 1-dof 
%   models),NgdlTot degrees of freedom, proportional coefficients alfa
%   (inertia) and beta (stiffness) and stiffness k1 for the first segment
%   of drill-pipes.
%
%   Inputs:
%   I           -> Inertia matrix
%   K           -> Stiffness matrix
%   Xi          -> Damping coefficient
%   NgdlTot     -> Total number of degrees of freedom at the drill string
%   alfa        -> Damping constant proportional to inertia
%   beta        -> Damping constant proportional to stiffness
%   k1          -> Stiffness of the first pipe segment
%   i_print     -> If True, print usefull information
%
%   Outputs:
%   C           -> Damping matrix
%   c1          -> Damping of the first pipe segment
%
%   Ex:
%   alfa = 0.5, beta = 0.5, Xi = 1, NgdlTot = 3
%   K = [200  -100   0;
%        -100 200 -100;
%         0   -100 100;
%   I = [200  0    0;
%         0  200   0;
%         0   0  200;
%   [C,c1] = C_Gen(I,K,Xi,NgdlTot,alfa,beta) 
%       
%   C = [200  -50   0;
%        -50  200 -50;
%         0   -50 150;
%
%   c1 = 50;
%
%   IMPORTANT!! To consider more dofs at the BHA, it is recommended to
%   increase the dofs at the DP as well.
%
%   LAST MODIFIED: 30/03/2020 BY JORDAN BARBOZA AND DANIEL LOBO
%   CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


[Eigen_vec,Eigen_val] = eig(K,I);
if length(I) == 1
    C = alfa*I + beta*K;
    Xi = 0.5 * Eigen_vec' * C * Eigen_vec * (sqrt(Eigen_val))^(-1);
    c1 = C;
else
    C = alfa*I + beta*K;
    Xi = 0.5 * Eigen_vec' * C * Eigen_vec * (sqrt(Eigen_val))^(-1);
    c1 = beta*k1;
end

if i_print
    disp('Torsional damping Matrix'),disp(C)
    disp('Torsional damping coefficients'),disp(diag(Xi))
    disp('Torsional natural frequencies (Hz):'),disp(sqrt(...
                                                     diag(Eigen_val))/2/pi)
end
end