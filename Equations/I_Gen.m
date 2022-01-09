function I = I_Gen(NgdlTotDP,NgdlTotBHA,IDP,IBHA,IBIT)

% I_Gen Generates the inertia tensor for a torsional lumped parameter
%       model of a drill-string considering the BHA as a rigid or flexible
%       body.
%
%   I = I_Gen(NgdlTotDP,NgdlTotBHA,IDP,IBHA,IBIT)
%   with I as the inertia matrix generated for a model with NgdlTotDP
%   degrees of freedom at the drill-pipes and NgdlTotBHA degrees of freedom
%   at the BHA, drill-pipes with moment of inertia IDP, BHA with moment of
%   inertia IBHA, drill-bit with moment of inertia IBIT.
%
%   Inputs:
%   NgdlTotDP   -> Total number of degrees of freedom at the drill pipes
%   NgdlTotBHA  -> Total number of degrees of freedom at the drill collars
%                  If NgdlTotBHA == 0, BHA is a rigid body
%   IDP         -> Drill pipes moment of inetia
%   IBHA        -> Drill collars moment of inertia 
%   IBIT        -> Bit bounce moment of inertia
%
%   Outputs:
%   I           -> Inertia matrix 
%
%   Ex:
%   NgdlTotDP = 5, NgdlTotBHA = 0, IDP = [I1 I2 I3 I4 I5], IBHA = [I6 I7]
%   I = I_Gen(NgdlTotDP,NgdlTotBHA,IDP,IBHA,IBIT) 
%   I = [(I1+I2)/2     0          0         0                0            ;
%            0     (I2+I3)/2      0         0                0            ;
%            0         0     (I3+I4)/2      0                0            ;
%            0         0          0    (I4+I5)/2             0            ;
%            0         0          0         0     I5/2 + sum(IBHA) + IBIT ]      
%   
%      
%   IMPORTANT!! To consider more dofs at the BHA, it is recommended to
%   increase the dofs at the DP as well.
%
%   LAST MODIFIED: 30/03/2020 BY JORDAN BARBOZA AND DANIEL LOBO
%   CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


% vector with the moments of inertia of each segment of the drill string
Ivector = vertcat(IDP,IBHA); 

if NgdlTotBHA > 0       %(Flexible BHA)
    Ngdl = NgdlTotDP + NgdlTotBHA;
else % NgdlTotBHA = 0   %(Rigid BHA)
    Ngdl = NgdlTotDP;
end

I = zeros(Ngdl,Ngdl);
for i = 1:Ngdl-1
    I(i,i) = (1/2*Ivector(i,1) + 1/2*Ivector(i+1,1));
end

if NgdlTotBHA == 0         %(Rigid BHA)
    if NgdlTotDP == 1      % Equivalent kinetic energy method
        I(Ngdl,Ngdl) = 1/3*Ivector(Ngdl,1) + sum(IBHA) + IBIT;
    else %NgdlTotDP > 1    % Proportional inertia method
        I(Ngdl,Ngdl) = 1/2*Ivector(Ngdl,1) + sum(IBHA) + IBIT;
    end
else %NgdlTotBHA > 0       %(Flexible BHA)
    I(Ngdl,Ngdl) = 1/2*Ivector(Ngdl,1) + IBIT;
end

disp('Inertia matrix'), disp(I)
end