function [dB_dx_Kx] = systemJacobianStateFeedback(x, K)

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
K1 = mb + 2*mw + 2*Iw1/(r*r);   
K2 = -Ib3 + Ib1 + mb*l*l;
K3 = Ib2 + mb*l*l;
K4 = mb*l;
theta=x(3);
rho1 = K4*K4*sin(theta)*sin(theta) + mb*Ib2 + 2.*(mw + Iw1/(r*r))* K3;
rho2 = Ib3 + 2*Iw2 + 0.5*d*d*(mw + Iw1/(r*r)) + K2*sin(theta)*sin(theta);
drho1_dtheta = 2*K4*K4*sin(theta)*cos(theta);
drho2_dtheta = 2*K2*sin(theta)*cos(theta);

%% 1) calculate partial derivatives of input-matrix B
dB2_dtheta = (-rho1*K4*sin(theta) - (K3/r + K4*cos(theta))*drho1_dtheta)/(rho1^2);
dB4_dtheta = (rho1*K4*sin(theta)/r + (K1 + K4*cos(theta)/r)*drho1_dtheta)/(rho1^2);
dB6_dtheta = (2*d*r*drho2_dtheta)/(rho2^2);

%% 1) define Jacobian entries
dB_dx_Kx = zeros(n,n);

dB_dx_Kx(2,3) = dB2_dtheta*(x(1)*(K(1,1)+K(2,1)) + x(2)*(K(1,2)+K(2,2)) + ...
                            x(3)*(K(1,3)+K(2,3)) + x(4)*(K(1,4)+K(2,4)) + ...
                            x(5)*(K(1,5)+K(2,5)) + x(6)*(K(1,6)+K(2,6)) );
                        
dB_dx_Kx(4,3) = dB4_dtheta*(x(1)*(K(1,1)+K(2,1)) + x(2)*(K(12)+K(2,2)) + ...
                            x(3)*(K(1,3)+K(2,3)) + x(4)*(K(1,4)+K(2,4)) + ...
                            x(5)*(K(1,5)+K(2,5)) + x(6)*(K(1,6)+K(2,6)) );
                        
dB_dx_Kx(6,3) = dB6_dtheta*(x(1)*(K(1,1)-K(2,1)) + x(2)*(K(1,2)-K(2,2)) + ...
                            x(3)*(K(1,3)-K(2,3)) + x(4)*(K(1,4)-K(2,4)) + ...
                            x(5)*(K(1,5)-K(2,5)) + x(6)*(K(1,6)-K(2,6)) );
end