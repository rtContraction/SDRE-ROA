function [A,B] = systemSDC(z)
%% 1) define WIP parameters:
g    =   9.81;            %gravitational acceleration in m/s
mb   = 277.000*1e-3;      %mass of the pendulum body in kg
mw   =  28.000*1e-3;      %mass of an individual wheel in kg
Ib1  = 543.108*1e-6;      %moment of inertia (moving direction) of the pendulum body in kg/m^2 
Ib2  = 481.457*1e-6;      %moment of inertia (in moving plane perpendicular) of the pendulum body in kg/m^2 
Ib3  = 153.951*1e-6;      %moment of inertia (out of moving plane perpendicular) of the pendulum body in kg/m^2 
Iw1  =   7.411*1e-6;      %moment of inertia of the wheel around its spinning axis in kg/m^2 
Iw2  =   4.957*1e-6;      %moment of inertia of the wheel perpendicular to spinning axis in kg/m^2 
l    =  48.670*1e-3;      %height of centre of mass of the body in upright position in m
r    =  33.100*1e-3;      %wheel radius in m
d    =  98.000*1e-3;      %distance between two wheels in m

%% 2) state variables initialisation:
%x = z(1);        %straight forward distance of wheeled inverted pendulum
xdot = z(2);      %straight forward velocity of wheeled inverted pendulum
theta = z(3);     %pitch angle
thetadot = z(4);  %angular velocity of pitch
%psi = z(5);      %yaw angle 
psidot = z(6);    %angular velocity of yaw

%% 3) introduce constants for more compact formulation of matrix entries
K1 = mb + 2*mw + 2*Iw1/(r*r);
K2 = -Ib3 + Ib1 + mb*l*l;
K3 = Ib2 + mb*l*l;
K4 = mb*l;

nu1 = K4*K4*sin(theta)*sin(theta) + mb*Ib2 + 2.*(mw + Iw1/(r*r))* K3;
nu2 = Ib3 + 2*Iw2 + 0.5*d*d*(mw + Iw1/(r*r)) + K2*sin(theta)*sin(theta);  

%% 4) set up B(z) matrix
B = zeros(6,2);
B(2,1) = (K3/r + K4*cos(theta)) * (1./nu1);
B(2,2) = B(2,1);
B(4,1) = ((K4*cos(theta))/r + K1) * (-1./nu1);
B(4,2) = B(4,1);
B(6,1) = (0.5*d/r) * (-1.0/nu2);                                          
B(6,2) = -B(6,1); 

%% 5) define A(z) matrix
%% 5.1) if theta ~ 0.0 for numerical stability sin(theta)/(theta) has to be
%% set to sin(theta)/(theta) => 1.0
%if (theta < 1e-14 && theta > -1e-14)
if (theta==0)
    theta_nu1 = 1.0/nu1;
else
    theta_nu1 = sin(theta)/(theta*nu1);
end

%% 5.2) compute matrix entries
A = zeros(6,6);
A(1,2) = 1;
A(2,3) = theta_nu1 * -(K4 * K4 * g * cos(theta));
A(2,4) = (sin(theta)/nu1) * K4 * K3 * thetadot;
A(2,6) = (sin(theta)/nu1) * K4 * (K3 - K2 * cos(theta) * cos(theta))*psidot;
A(3,4) = 1;
A(4,3) = theta_nu1 * (K4 * K1 * g);
A(4,4) = (sin(theta)/nu1) * -(K4 * K4 * cos(theta)) * thetadot;
A(4,6) = (sin(theta)/nu1) * (-K4 * K4 + K1 * K2) * cos(theta) * psidot;
A(5,6) = 1;
A(6,2) = (sin(theta)/nu2) * -(0.5 * K4) * psidot;
A(6,4) = (sin(theta)/nu2) * -(K2 * cos(theta)) * psidot;
A(6,6) = (sin(theta)/nu2) * -(0.5 * K4 * xdot + K2 * cos(theta) * thetadot);
end