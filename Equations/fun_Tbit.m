function Tbit = fun_Tbit(vphi,WOB)
%  fun_Tbit  Calculates the torque on bit according to a pure torsional 
%            regularized bit-rock interaction model
% 
%   Tbit = fun_Tbit(vphi,WOB)
%   with Tbit as the torque on bit
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
%     WOBref = 244.2e+3;
%     muref = 1.4;
%     b0 = 4.17e+3;
%     b1 = 1.91;
%     b2 = 8.5;
%     b3 = 5.47;
%     Xt = 0.0239

Tbit =  WOB*0.0239*(tanh(1.91*vphi) + 8.5*vphi./(1+5.47*vphi.^2));
       
end