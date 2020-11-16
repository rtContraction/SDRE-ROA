function [dA_dx_x] = systemJacobianOpenLoop(x)

%% 0) Define function handles for Nonlinear Terms
n    =   6;
g    =   9.81;                  %gravitational acceleration in m/s
mb   = 277.000*1e-3;            %mass of the pendulum body in kg
mw   =  28.000*1e-3;            %mass of an individual wheel in kg
Ib1  = 543.108*1e-6;            %moment of inertia (moving direction) of the pendulum body in kg/m^2 
Ib2  = 481.457*1e-6;            %moment of inertia (in moving plane perpendicular) of the pendulum body in kg/m^2 
Ib3  = 153.951*1e-6;            %moment of inertia (out of moving plane perpendicular) of the pendulum body in kg/m^2 
Iw1  =   7.411*1e-6;            %moment of inertia of the wheel around its spinning axis in kg/m^2 
Iw2  =   4.957*1e-6;            %moment of inertia of the wheel perpendicular to spinning axis in kg/m^2 
l    =  48.670*1e-3;            %height of centre of mass of the body in upright position in m
r    =  33.100*1e-3;            %wheel radius in m
d    =  98.000*1e-3;            %distance between two wheels in m

%% 1.1) initialise states
xdot = x(2);
theta = x(3);
thetadot = x(4);
psidot = x(6);

%% 1.2) introduce constants for easier computation
K1 = mb + 2*mw + 2*Iw1/(r*r);   
K2 = -Ib3 + Ib1 + mb*l*l;
K3 = Ib2 + mb*l*l;
K4 = mb*l;

rho1 = K4*K4*sin(theta)*sin(theta) + mb*Ib2 + 2.*(mw + Iw1/(r*r))* K3;
rho2 = Ib3 + 2*Iw2 + 0.5*d*d*(mw + Iw1/(r*r)) + K2*sin(theta)*sin(theta);  
drho1_dtheta = 2*K4*K4*sin(theta)*cos(theta);
drho2_dtheta = 2*K2*sin(theta)*cos(theta);

%% 1) compute derivatives of A-entries
%% 1.1) theta derivatives
if (theta==0)
    dA43_dtheta = 0;
    dA23_dtheta = K4*K4*g/(rho1^2);  % will vanish
else
dA23_dtheta = (K4*K4*g/(rho1*theta)) * ( sin(theta)*sin(theta) - ...
               cos(theta)*cos(theta) + sin(theta)*cos(theta)*...
               (drho1_dtheta*theta + rho1)/(rho1*theta) );
dA43_dtheta = K1*K4*g*(rho1*(theta*cos(theta)-sin(theta))...
    -  theta*drho1_dtheta*sin(theta))...
              /(rho1*rho1*theta*theta); 
end

dA24_dtheta = K4*K3*thetadot*(rho1*cos(theta) - sin(theta)*drho1_dtheta)...
              /(rho1*rho1);
dA26_dtheta = K4*psidot*(rho1*(K3*cos(theta) - K2*cos(theta)^3 + ...
              2*K2*cos(theta)*sin(theta)^2) - sin(theta)*drho1_dtheta*...
              (K3 - K2*cos(theta)^2)) /(rho1*rho1);      
dA44_dtheta = K4*K4*thetadot*((sin(theta)^2 - cos(theta)^2)*rho1 + ...
              sin(theta)*cos(theta)*drho1_dtheta)/(rho1*rho1);     
dA46_dtheta = (K1*K2 - K4*K4)*psidot*((- sin(theta)^2 + cos(theta)^2)*rho1 - ...
              sin(theta)*cos(theta)*drho1_dtheta)/(rho1*rho1);  
dA62_dtheta = 0.5*K4*psidot*(sin(theta)*drho2_dtheta - rho2*cos(theta))...
              /(rho2*rho2);
dA64_dtheta = K2*psidot*(sin(theta)*cos(theta)*drho2_dtheta + rho2*...
              (sin(theta)^2 - cos(theta)^2) )/(rho2*rho2);
dA66_dtheta = 0.5*K4*xdot*(sin(theta)*drho2_dtheta - rho2*cos(theta))...
              /(rho2*rho2) + ...
              K2*thetadot*(sin(theta)*cos(theta)*drho2_dtheta + rho2*...
              (sin(theta)^2 - cos(theta)^2) )/(rho2*rho2);
          
%% 1.2) thetadot derivatives
dA24_dthetadot = (sin(theta)/rho1)*K3*K4;          
dA44_dthetadot = (sin(theta)/rho1)*-K4*K4*cos(theta); 
dA66_dthetadot = (sin(theta)/rho2)*-K2*cos(theta); 

%% 1.3) psidot derivatives
dA26_dpsidot = (sin(theta)/rho1)*K4*(K3 - K2*cos(theta)^2);          
dA46_dpsidot = (sin(theta)/rho1)*(K1*K2 - K4*K4)*cos(theta); 
dA62_dpsidot = (sin(theta)/rho2)*-0.5*K4;
dA64_dpsidot = (sin(theta)/rho2)*-K2*cos(theta);

%% 1.4) xdot derivatives
dA66_dxdot = (sin(theta)/rho2)*-0.5*K4;

%% 2) initialise Jacobian and set non-zero-entries
dA_dx_x = zeros(n,n);

dA_dx_x(2,3) = dA23_dtheta*theta + dA24_dtheta*thetadot + dA26_dtheta*psidot;
dA_dx_x(2,4) = dA24_dthetadot*thetadot;
dA_dx_x(2,6) = dA26_dpsidot*psidot;

dA_dx_x(4,3) = dA43_dtheta*theta + dA44_dtheta*thetadot + dA46_dtheta*psidot;
dA_dx_x(4,4) = dA44_dthetadot*thetadot;
dA_dx_x(4,6) = dA46_dpsidot*psidot;

dA_dx_x(6,2) = dA66_dxdot * psidot;
dA_dx_x(6,3) = dA62_dtheta*xdot + dA64_dtheta*thetadot + dA66_dtheta*psidot;
dA_dx_x(6,4) = dA66_dthetadot*psidot;
dA_dx_x(6,6) = dA62_dpsidot*xdot + dA64_dpsidot*thetadot;
end