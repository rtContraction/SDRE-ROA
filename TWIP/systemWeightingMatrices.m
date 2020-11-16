function [Q, R] = systemWeightingMatrices(z)
%% given a state-vector z, this script computes the according weighting 
%% matrices Q(z) and R(z) to define the performance index for the 
%% SDRE-controller
%% 1) state variables initialisation:
x = z(1);        %straight forward distance of wheeled inverted pendulum
xdot = z(2);     %straight forward velocity of wheeled inverted pendulum
theta = z(3);    %pitch angle
thetadot = z(4); %angular velocity of pitch
%psi = z(5);      %yaw angle 
psidot = z(6);   %angular velocity of yaw

%% 2) define R(z) matrix for input penalization 
%% phyiscal motivation: psidto is only really critical 
%R_const = 0.5.*eye(2);

% R = eye(2).*(0.5 + 400.0*(psidot^4) + 150.0*thetadot*thetadot);

fz=40*psidot^2 + 40*thetadot*thetadot+ 60*xdot*xdot;
%fz0=4*(4*pi)^2+4*(4*pi)^2+4*(2^2);
fz0=2*(4*pi)^2+2*(4*pi)^2+2*(2^2);
g=(exp(-fz/fz0));
R = eye(2).*(1e5-9.95e4*g);
 %R = eye(2).*(1e5);
%% 3) define Q(z) matrix for state penalization 
Q = zeros(6,6);
% Q_const(1,1) = 1;                                                             
% Q_const(2,2) = 3.0;                          
% Q_const(3,3) = 1.0;
% Q_const(4,4) = 3.0;
% Q_const(5,5) = 1.0;
% Q_const(6,6) = 2.5; 
% 
% Q=Q_const;
% 
% Q(1,1) = 4.5;                                                             
% Q(2,2) = 3.0  + 400.0*(theta*theta/pi*pi)*psidot*psidot +... 
%                 200.0*(theta*theta/pi*pi)*xdot*xdot;                   
% Q(3,3) = 10.0 + 2000.0*(theta*theta/pi*pi);
% Q(4,4) = 2.0  + 1500.0*(theta*theta/pi*pi)*(thetadot*thetadot) +...
%                 500.0*(theta*theta/pi*pi)*(psidot*psidot);
% Q(5,5) = 4.0;
% Q(6,6) = 1.0  + 100.0*(theta*theta/pi*pi)*(thetadot*thetadot + xdot*xdot) + ...
%                 20.0*(theta*theta/pi*pi)*psidot*psidot;
% k_gain=.5e0;
% Q(1,1) = 100 + k_gain*200.0*exp(-xdot*xdot) ;  
% 
% Q(2,2) = 3.0  + k_gain*2.0*(theta*theta/pi*pi)*psidot*psidot +... 
%                  k_gain*1.0*(theta*theta/pi*pi)*xdot*xdot;                   
% Q(3,3) = 1.0 +  k_gain*4.0*(theta*theta/pi*pi);
% Q(4,4) = 2.0  +  k_gain*1.5*(theta*theta/pi*pi)*(thetadot*thetadot) +...
%                  k_gain*5.0*(theta*theta/pi*pi)*(psidot*psidot);
% Q(5,5) = 100;
% Q(6,6) = 1.0  +  k_gain*1.0*(theta*theta/pi*pi)*(thetadot*thetadot + xdot*xdot) + ...
%                  k_gain*.20*(theta*theta/pi*pi)*psidot*psidot;

k_gain=.5e0;
Q(1,1) = 100 + k_gain*200.0*exp(-xdot*xdot) ;  

Q(2,2) = 3.0  + k_gain*2.0*(theta*theta/pi*pi)*psidot*psidot +... 
                 k_gain*1.0*(theta*theta/pi*pi)*xdot*xdot;                   
Q(3,3) = 1.0 + k_gain*200.0*exp(-thetadot*thetadot);
Q(4,4) = 2.0  +  k_gain*1.5*(theta*theta/pi*pi)*(thetadot*thetadot) +...
                 k_gain*5.0*(theta*theta/pi*pi)*(psidot*psidot);
Q(5,5) = 100 + k_gain*200.0*exp(-psidot*psidot);
Q(6,6) = 1.0  +  k_gain*1.0*(theta*theta/pi*pi)*(thetadot*thetadot + xdot*xdot) + ...
                 k_gain*.20*(theta*theta/pi*pi)*psidot*psidot;             
             

end
