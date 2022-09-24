function dydt = ode_4gdl(t,y,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
                     kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall)
% ode_4gdl  Drill-string system of equation. Calculates the time derivative
%           of the state space for a torsional lumped parameter model of a
%           drill-string.
% 
%   dydt = ode_4gdl(t,y,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,kt,ct,...
%                                                             Omega,WOB)
%   with dydt as the time derivative of the state space.
% 
%   Inputs:
%   t         -> Initial time
%   y         -> Initial conditions
%   e         -> Eccentricity
%   u         -> Wall friction coefficient
%   ks        -> Borehole wall stiffness
%   cs        -> Borehole wall damping
%   R         -> Collar Radius
%   ch        -> Fluid drag damping
%   k         -> Stiffness (N/m)
%   tol       -> Gap between collar and borehole wall
%   Lc        -> BHA section length
%   Mt        -> Total mass (kg)
%   aphi      -> Initial acceleration
%   Im        -> Inertia matrix
%   kt        -> Torsional stiffness of each ktp element
%   ct        -> Proportional Damping
%   Omega     -> Rotational speed in rad/s
%   WOB       -> Weight on Bit in [N]
%
%   Outputs:
%   dydt      -> Time derivative of the state space
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


% Extract system's size
lp_number = length(Lc);

% Allocate memory 
H     = zeros(lp_number,1); % impact boolean
Fn    = zeros(lp_number,1); % impact normal force
Ft    = zeros(lp_number,1); % impact friction force
Tlati = zeros(lp_number,1); % i-th section's impact torque
dydt  = zeros(length(y),1); % State-space vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Torque on bit
Tbit = fun_Tbit(y(2*length(Im)),WOB); 

i = 2*length(Im)+1;
for j = 1:lp_number

    % Stochastic field
    H_s = Hs_extract(H_grid,z_grid,theta_grid,y(i+2),y(end));
    tol(j) = tol(j) + H_s*Dbwall/2;
    
    %Check Impact
    if y(i) > tol(j)
        H(j) = 1;
    else
        H(j)=0;
    end
    
        % Impact normal force
    Fn(j)    = fun_Fr(y(i),y(i+1),cs,ks,H(j),tol(j)); 
    
    % Impact friction
    Ft(j)    = fun_Ftheta(y(2*length(Im)),y(i),y(i+3),u,Fn(j),R(j)); 
    
    % j-th section's impact torque
    Tlati(j) = fun_Tlat(y(length(Im)),y(i),y(i+1),y(i+2),y(i+3),...
        ch(j),e(j),R(j),Fn(j),Ft(j));   
    i = i + 4;
end

% Total impact torque
Tlat         = sum(Tlati); 

%%%%%%%%%%%%%%%%%%%% CALCULATE STATE SPACE DERIVATIVE %%%%%%%%%%%%%%%%%%%%%

dydt(1:2*length(Im),1)  = torcional(t,y(1:length(Im)),...
    y(length(Im)+1:2*length(Im)),kt,ct,k1,c1,Omega,Tbit,Im,invIm,Tlat);

i = 2*length(Im)+1;
for j = 1:lp_number 
    dydt(i:i+3,1)   = lateral(y(length(Im)),y(2*length(Im)),y(i),...
        y(i+1),y(i+2),y(i+3) ,aphi,...
        k(j),Mt(j),Lc(j),Tbit,e(j),ch(j),Ft(j),Fn(j));
    i = i + 4;
end

dydt(end,1) = fun_ROP(y(2*length(Im)),WOB);
 

end