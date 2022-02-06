function [V,w,MN] = FreqNat_Gen(K,I,NgdlTotDP,NgdlTotBHA,i_norm)

% FreqNat_Gen Generates natural frequencies, eigenvectors and normal modes
%             for a torsional lumped parameter model of a drill-string
%             considering the BHA as a rigid or flexible body.
%
%   [V,w,MN] = FreqNat_Gen(K,I,NgdlTotDP,NgdlTotBHA,NgdlTot,i_norm) 
%   returns the solution of the problem K*V = w^2*M*V with V, w and MN as 
%   the eigenvectors, natural frequencies and normal modes, respectively,
%   generated for a model with NgdlTotDP degrees of freedom at the 
%   drill-pipes and NgdlTotBHA degrees of freedom at the BHA
%
%   Inputs: 
%   K          -> Stiffness matrix
%   I          -> Inertia matrix
%   NgdlTotDP  -> Total number of degrees of freedom at the drill pipes
%   NgdlTotBHA -> Total number of degrees of freedom at the drill collars
%                 If NgdlTotBHA == 0, BHA is a rigid body
%   i_norm     -> Determines the normalization of plotting of the modes
%        
%   Outputs:
%   w          -> Natural frequencies
%   V          -> Eigenvectors
%   MN         -> Normal modes
%   
%   Ex:
%   K =  1.0e+03 *[1.9485   -0.9742         0;
%                 -0.9742    1.9485   -0.9742;
%                     0     -0.9742    0.9742];
%   I = [221.8703         0         0;
%          0           221.8703     0;
%          0              0      344.6975];
%   NgdlTotDP = 3; NgdlTotBHA = 0; i_norm = 1
%   [V,w,MN] = FreqNat_Gen(K,I,NgdlTotDP,NgdlTotBHA,NgdlTot,i_norm)
%
% V = [0.0183   -0.0481   -0.0431;
%      0.0338   -0.0310    0.0490;
%      0.0442    0.0281   -0.0126]
% 
% w = [0.1298;
%      0.3882;
%      0.5907]
% 
% MN = [0.4137    1.0000   -0.8792;
%       0.7647    0.6451    1.0000;
%       1.0000   -0.5838   -0.2581]
%    
%   LAST MODIFIED: 10/03/2020 BY JORDAN BARBOZA AND DANIEL LOBO
%   CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


[V,D] = eig(K,I);
EigenValues = sqrt(D)/2/pi;
FreqNatHzVetor = diag(EigenValues);
w = FreqNatHzVetor;

if NgdlTotBHA > 0       %(Flexible BHA)
    Ngdl = NgdlTotDP + NgdlTotBHA;
else % NgdlTotBHA = 0   %(Rigid BHA)
    Ngdl = NgdlTotDP;
end

ModosNormais = zeros(Ngdl,Ngdl);
for i = 1:Ngdl
    max_min = [max(V(:,i)) min(V(:,i))];
    if i_norm == 1 % Set the higher amplitude to 1
        [fator,kk] = max(abs(max_min));
        fator = fator*sign(max_min(kk));
    end
    if i_norm == 0 % Set the drill amplitude to 1
        fator = V(end,i); 
    end
    ModosNormais(:,i) = V(:,i)/fator;
end
MN = ModosNormais;
end