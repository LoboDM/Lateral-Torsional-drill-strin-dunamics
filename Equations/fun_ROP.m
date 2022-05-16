function ROP = fun_ROP(vphi,WOB)
%  fun_ROP  Calculates the rate of penetration using a regularized 
%           bit-rock interaction model
% 
%   ROP = fun_ROP(vphi,WOB)
%   with ROP as the rate of penetration
% 
%   Inputs:
%   vphi -> Drill-bit speed
%   WOB  -> Weight on Bit in [N]
%
%   Outputs:
%   Tbit -> Torque on bit
% 
%  LAST MODIFIED: 11/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


% Reference constants:       
%     a1 = 3.429e-3;
%     a2 = 5.672e-8;
%     a3 = 1.374e-4;

a1 = 3.429e-3;
a2 = 5.672e-8;
a3 = 1.374e-4;

Z = vphi./sqrt(vphi.^2 + 2^2);
ROP = Z.^2*(-a1 + a2*WOB) + Z.*a3.*vphi;
       
end