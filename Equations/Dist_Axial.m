function Fa = Dist_Axial(WOB,Lp,Dpo,Dpi,Lc,Dco,Dci,rho,g)
% Dist_Axial  Calculate axial force in the middle of each lumped parameter.
% 
%   Fa = Dist_Axial(WOB,Lp,Dpo,Dpi,Lc,Dco,Dci,rho,g)
%   with Fa as the axial force in each section.
% 
%   Inputs:
%   WOB -> Weight on Bit in [N]
%   Lp  -> Drill pipes length
%   Dpo -> Drill pipes outer diameter
%   Dpi -> Drill pipes inner diameter
%   Lc  -> BHA section length
%   Dco -> BHA collar outer diameter
%   Dci -> BHA collar inner diameter
%   rho -> Specific mass
%   g   -> Acceleration of gravity
%
%   Outputs:
%   Fa  -> Axial force in each section
% 
%  LAST MODIFIED: 08/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


Lc  = [Lp Lc];
Dco = [Dpo Dco];
Dci = [Dpi Dci];
lLc = length(Lc);
Lt  = sum(Lc);

A = zeros(lLc,1);
dWOB = zeros(lLc,1);
Fa   = zeros(lLc,1);

for i = 1:lLc
    A(i) = pi*(Dco(i)^2-Dci(i)^2)/4;
    dWOB(i) = rho*A(i)*g;
end

% Creates force in the form of Fa = dWOB x + b, where x is the length.
b = zeros(lLc,1);
L = Lt;
b(end) = WOB - dWOB(end)*L;
for i = (lLc-1):-1:1
    L = Lt - Lc(end-i);
    b(i) = (dWOB(i+1)-dWOB(i))*L + b(i+1);
end

% Calculates axial force in the middle of the section
for i = 1:lLc
    x     = Lt - Lc(i)/2 - sum(Lc(i+1:end)); 
    Fa(i) = dWOB(i)*x + b(i);
end

% Remove axial force from the drill-pipes
Fa(1) = [];
