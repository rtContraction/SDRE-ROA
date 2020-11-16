function [sdre,sys]=systemConfig(box,points,degree,Const)
SDRE_Restart % clean variables, including YALMIP
% ==================  calculates the Metric M(x)===========================
[M,E,Md,trajBound,sys] = systemData(box,points,degree,Const);
% ===============  Default Parameters for  the Simulation =================
sys.cRateChang=20;
sys.lambda_min=trajBound(2);
sys.lambda_max=trajBound(3);
sys.overS=sys.lambda_max/sys.lambda_min;
sys.box=box;
sys.Const=Const;
sys.B=[0;1];
sys.A=@(x) [x(1) 1;0 0];
sys.J=@(x) [2*x(1)  1;0 0];
% ===============   SDRE structure default values  ========================
sdre.n=2;        % states
sdre.m=1;        % inputs
sdre.x0=[0;0];   % Initial condition for the system
sdre.Q=100*eye(2);
sdre.R=1;
sdre.comp=2;                   
sdre.simT=5;     % simulation time
sdre.ts=1e-2;    % simulation step
sdre.M=M;        % Metric from  optimisation
sdre.Md=Md;      % dM/dt
sdre.E=E;        % symbolic representation of  norm(x)_M^2 = E = xT *M(x)*x
sdre.cRate=trajBound(1);
% ===================  Transformation of Coordinates  =====================
sys.Tx = @(x) x(:);
sys.Tu = @(u) u(:);
% ===================  AntiTransformation of Coordinates  =================
sys.aTx=@(x) x(:);
sys.aTu=@(u) u(:);
% ===================  Stationary State  ==================================
sys.X_ss=[0;0];
sys.U_ss=0;
% ===============   Function to calculate ur ==============================
sdre.track=0;
end




