function fun_Model(ti,tf,dt,tolerance,rho,rho_f,Dco,Dci,Ca,Dpo,Dpi,...
                Dbwall,Lp,E,G,Cd,g,u,ks,cs,alpha,a_c,b_c,Lc,rpm,WOBf,...
                local,LATERAL_dofs,N_tor,i_print)
% fun_Model   Generates the analysis and simulation of a lateral-torsional
%             lumped parameter model of a drill-string considering the
%             axial force in stiffness matrices.  
%
%   Inputs:
%   ti        -> Initial time                                   
%   tf        -> Final time
%   dt        -> Maximum step/Saving step                       
%   tolerance -> Solver tolerance                               
%   rho       -> Specific mass                                  
%   rho_f     -> Fluid specific mass                            
%   Dco       -> BHA collar outer diameter
%   Dci       -> BHA collar inner diameter
%   Ca        -> Added mass coefficient
%   Dpo       -> Drill pipes outer diameter
%   Dpi       -> Drill pipes inner diameter
%   Dbwall    -> Borehole wall diameter
%   Lp        -> Drill pipes length
%   E         -> Modulus of elasticity
%   G         -> Shear coefficient
%   Cd        -> Drag coefficient
%   g         -> Acceleration of gravity
%   u         -> Wall friction coefficient
%   ks        -> Borehole wall stiffness
%   cs        -> Borehole wall damping
%   alpha     -> Drill-string slope
%   a_c       -> Damping constant proportional to inertia
%   b_c       -> Damping constant proportional to stiffness
%   Lc        -> BHA section length
%   rpm       -> Rotational speed of rotary table in rpm
%   WOBf      -> Axial force at the bit
%   local     -> Save folder
%   LATERAL_dofs-> Indices from Lc considered in lat. dyn.
%   N_tor      -> Number of torsional DOfs.           
%
%  Outputs:
%  analysis and simulation of a torsional lumped parameter model
% 
%  LAST MODIFIED: 11/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


%%%%%%%%%%%%%%%%%%%%%%%% LUMPED PARAMETERS TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotational speed in rad/s
Omega = rpm*2*pi/60;

% Axial forces
[WOB]= Dist_Axial(WOBf,Lp,Dpo,Dpi,Lc,Dco,Dci,rho,g);
lLc = length(Lc);  
for i = lLc:-1:1
    
    % Lp mass
    m(i)       = pi*rho*(Dco(i)^2-Dci(i)^2)*Lc(i)/8;        
    
    % Added fluid mass
    mf(i)      = pi*rho_f*(Dci(i)^2 + Ca*Dco(i)^2)*Lc(i)/8; 
    
    % Total mass (kg)
    Mt(i)      = m(i)+mf(i);                               
    
    %Area intertia (m4)
    I_area(i)  = pi*(Dco(i)^4-Dci(i)^4)/64;                 
    
    % Stiffness (N/m)
    k(i)       = (pi^4)*E*I_area(i)/(2*Lc(i)^3) - WOB(i)*pi^2/(2*Lc(i)); 
    
    % Fluid drag damping
    ch(i)      = 2*(rho_f*Cd*Dco(i)*Lc(i)/(3*pi)); 
    
    % Collar Radius
    R(i)       = Dco(i)/2;                  
    
    % Gap between collar and borehole wall
    tol(i)     = (Dbwall -Dco(i))/2;        
    
    % Eccentricity
    e(i)       = Dpo/10;                    
    
    if k(i) < 0
        % Buckling checkage
        h = msgbox('k < 0 ' );
        return
    end
end


%%%%%%%%%%%%%%%%%%%% Calculated Values -- Torsional %%%%%%%%%%%%%%%%%%%%%%%
NgdlTotDP = N_tor*length(Lp);

% Drill pipes incidence matrix
inc_DP = zeros(NgdlTotDP,length(Lp));
jj = 1;
for ii = 1:N_tor:NgdlTotDP
    inc_DP(ii:ii+N_tor-1,jj) = 1;
    jj = jj + 1;
end

inc_BHA = eye(length(Lc));

DoutDPv = inc_DP*Dpo;
DinDPv = inc_DP*Dpi;
LDPv = inc_DP*Lp/N_tor;

DoutBHAv = inc_BHA*Dco';
DinBHAv = inc_BHA*Dci';
LBHAv = inc_BHA*Lc';

JDP = pi*(DoutDPv.^4-DinDPv.^4)/32;   % drill pipes polar moment of inertia
IDP = rho*JDP.*LDPv;     % drill pipes moment of inertia
kDP = G*JDP./LDPv;      % drill pipes stiffness

JBHA = pi*(DoutBHAv.^4-DinBHAv.^4)/32;% BHA polar moment of inertia
IBHA = rho*JBHA.*LBHAv;  % drill collars moment of inertia
kBHA = G*JBHA./LBHAv;   % drill collars stiffness

% Generates the stiffness matrix -> K
[kt,k1] = K_Gen(NgdlTotDP,0,kDP,kBHA,i_print);

% Generates the inertia matrix -> I
Im = I_Gen(NgdlTotDP,0,IDP,IBHA,0,i_print);
invIm = inv(Im); 

% Generates the damping matrix -> C
[ct,c1] = C_Gen(Im,kt,a_c,b_c,k1,i_print);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%
LATERAL_dofs = sort(LATERAL_dofs);
ci = zeros(4*length(LATERAL_dofs) + 2*length(Im),1);
jj = 1;
for i = LATERAL_dofs
    if m(i)*g*sin(alpha)/k(i) < tol(i)
        Ro(jj) = m(i)*g*sin(alpha)/k(i);
    else
        Ro(jj) = tol(i);
    end
    ci(2*length(Im) + 1 + (jj-1)*4) = Ro(jj);
    jj = jj + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Format and solve ode %%%%%%%%%%%%%%%%%%%%%%%%%%%%

aphi = 0;                   % Initial acceleration
tspan = [ti tf];            % Time information

[X,t] = rkf45drillstring(e(LATERAL_dofs),u,ks,cs,R(LATERAL_dofs),...
    ch(LATERAL_dofs),k(LATERAL_dofs),tol(LATERAL_dofs),Lc(LATERAL_dofs),...
    Mt(LATERAL_dofs),aphi,Im,invIm,kt,ct,k1,c1,Omega,WOBf,ci,dt,...
    tspan,tolerance);

% Save values
phi  = X(length(Im),:);
vphi = X(length(Im)*2,:);
lumped_parameter = length(LATERAL_dofs);
[lin,col] = size(phi);

r      = zeros(lumped_parameter,col);
vr     = zeros(lumped_parameter,col);
theta  = zeros(lumped_parameter,col);
vtheta = zeros(lumped_parameter,col);

j = 2*length(Im)+1;
for i = 1:lumped_parameter
    r(i,:)     = X(j  ,:);
    vr(i,:)    = X(j+1,:);
    theta(i,:)  = X(j+2,:);
    vtheta(i,:) = X(j+3,:);
    j = j + 4;
end

nome = strcat(local,'\WOB =',num2str(round(WOBf)),'rpm =',num2str(rpm,...
    '%03.f'));

save(nome);

end