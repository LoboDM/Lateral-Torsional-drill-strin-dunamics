function [X,tspan] = rkf45drillstring(e,u,ks,cs,R,ch,k,tol,Lc,Mt,...
    aphi,Im,invIm,kt,ct,k1,c1,Omega,WOB,xi,dt,tspan,tolerance,z_grid,...
    theta_grid,H_grid,Dbwall)
% rkf45drillstring  Runge-Kutta-Felbergh adaptative solver. It saves the
%                   acceleration vector to use as input in the system's
%                   nonlinear forces.
% 
%   [X,tspan] = rkf45drillstring(e,u,ks,cs,R,ch,k,tol,Lc,Mt,...
%                           aphi,Im,kt,ct,Omega,WOB,xi,dt,tspan,tolerance)
%   with X as the coordinate vector and tspan as the time information.
% 
%   Inputs:
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
%   xi        -> Initial conditions
%   dt        -> Maximum step/Saving step
%   tspan     -> Time information
%   tolerance -> Solver tolerance
%
%   Outputs:
%   X         -> Coordinate vector 
%   tspan     -> Time information
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


%%%%%%%%%%%%%%%%%%%%%%%%% METHOD CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RK5 Constants
c51 = 16/135;
c52 = 6656/12825;
c53 = 28561/56430;
c54 = -9/50;
c55 = 2/55;

% RK4 Constants
c41 = 25/216;
c42 = 1408/2656;
c43 = 2197/4104;
c44 = -1/5;

% Error Constants
ce1 = 1/360;
ce3 = -128/4275;
ce4 = -2197/75240;
ce5 = 1/50;
ce6 = 2/55;

% Constants for Ks
% K2
ck21 = 1/4; ck22 = 1/4;
% K3
ck31 = 3/8;  ck32 = 3/32; ck33 = 9/32;
% K4
ck41 = 12/13; ck42 = 1932/2197; ck43 = -7200/2197; ck44 = 7296/2197;
% K5
ck51 = 1; ck52 = 439/216; ck53 = -8; ck54 = 3680/513; ck55 = -845/4104;
% K6
ck61 = 1/2; ck62 = -8/27; ck63= 2; ck64 = -3544/2565; ck65 = 1859/4104;...
    ck66 = -11/40;


%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE SOLUTION VECTORS %%%%%%%%%%%%%%%%%%%%%%%%
ti = tspan(1);
tf = tspan(end);
h  = dt;
%
% t  = zeros(size(ti:dt:tf));
% X  = zeros(length(xi),length(t));
t(1) = ti;
[x_l,x_c] = size(xi);
if x_c > 1
    xi = xi';
end

tspan  = ti:h:tf;

X = zeros(length(xi),length(tspan));
X(:,1) = xi;
h_aux = h;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STARTS LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 2;
aphi_aux = aphi;
while tspan(i) <= tf
    ok = 0;
    atol  = tolerance;
    while ok == 0
        
        aphi_aux = aphi;
        K1 = h*ode_4gdl(ti,xi,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
                        kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
        aphi = K1(2*length(Im))/h;
        K2 = h*ode_4gdl(ti + ck21*h, xi + ck22*K1, e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
                        kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
        aphi = K2(2*length(Im))/h;
        K3 = h*ode_4gdl(ti + ck31*h, xi + ck32*K1 + ck33*K2,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
                        kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
        aphi = K3(2*length(Im))/h;
        K4 = h*ode_4gdl(ti + ck41*h, xi + ck42*K1 + ck43*K2 + ck44*K3,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
                        kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
        aphi = K4(2*length(Im))/h;
        K5 = h*ode_4gdl(ti + ck51*h, xi + ck52*K1 + ck53*K2 + ck54*K3 + ...
            ck55*K4,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
            kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
        aphi = K5(2*length(Im))/h;
        K6 = h*ode_4gdl(ti + ck61*h, xi + ck62*K1 + ck63*K2 + ck64*K3 + ...
            ck65*K4 - ck66*K5,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
            kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
        aphi = K6(2*length(Im))/h;
        
        % Fourth Order Runge-Kutta
        %         Xrk4  = xi  + (c41*K1  +  c42*K3  +  c43*K4  +  c44*K5);
        % Fifth Order Runge-Kutta
        Xrk5  = xi  + (c51*K1  +  c52*K3  +  c53*K4  +  c54*K5  +  c55*K6);
        % Estimated error
        error = norm(ce1*K1 + ce3*K3 + ce4*K4 + ce5*K5 + ce6*K6);
        % Tolerance for this step
        tolstep = tolerance*norm(xi) + atol;
        
        % Evaluates whether the error is less than the tolerance
        if error <= tolstep
            % Approves the time step
            ok     = 1;
            K1 = h*ode_4gdl(ti,xi,e,u,ks,cs,R,ch,k,tol,Lc,Mt,aphi,Im,invIm,...
                           kt,ct,k1,c1,Omega,WOB,z_grid,theta_grid,H_grid,Dbwall);
            aphi = K1(2*length(Im))/h;
            
            ti     = ti + h;
            xi     = Xrk5;
            
            if ti  == tspan(i)
                X(:,i) = Xrk5;
                i      = i + 1;
                if ti == tf, return, end
            end
            
        else
            aphi = aphi_aux;
            ok     = 0;
        end
        
        h      = h*(tolstep/error/2)^0.25;
        
        if ti + h > tspan(i)
            h      = tspan(i) - ti;
        elseif h > h_aux
            h      = h_aux;
        end
    end
end

end