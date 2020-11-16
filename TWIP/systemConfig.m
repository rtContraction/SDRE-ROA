 function [sdre,sys]=systemConfig(stateLimits, nPoints, polyGrad,Const,BisEnable)
SDRE_Restart % clean variables, including YALMIP

[Model, Mnum, E,Md, Coeffs, t, feasible, R_nonfeasible,...
 M_nonfeasible, MR_nonfeasible,sys] = systemData(...
stateLimits,nPoints,polyGrad,Const,BisEnable);
% ===============  CONTRACTION RATE PLOT SETTINGS =================
sys.cRate=t;
%sys.lambda_min=trajBound(2);
%sys.lambda_max=trajBound(3);
%sys.overS=sys.lambda_max/sys.lambda_min;
sys.box=stateLimits;
sys.feasible=feasible;
sys.Rnonfeasible=R_nonfeasible;
sys.Mnonfeasible=M_nonfeasible;
sys.MR=MR_nonfeasible;

% ===============   STATIC Simulation Parameters  ==================
sdre.n=6;        % states
sdre.m=2;        % inputs
sdre.comp=1;     % comp diff of 0 (this calls other type of controller
sdre.M=Mnum;     % Metric from  optimisation
sdre.Md=Md;      % dM/dt
sdre.E=E;        % symbolic representation of  norm(x)_M^2 = E = xT *M(x)*x
sdre.track=0;
% ===================  AntiTransformation of Coordinates  =================
% Path = [s,psi,theta,theta_dot,v,psi_dot]; u=[uR;uL];
% Kim =  [s,v,theta,theta_dot,psi,psi_dot]; u=[uL;uR];
sys.Tx = @(x) [x(1);x(5);x(3);x(4);x(2);x(6)]; % from Phatak -- > Kim
sys.Tu = @(u) [u(2);u(1)];
% ===================  AntiTransformation of Coordinates  =================
sys.aTx=@(x) [x(1);x(5);x(3);x(4);x(2);x(6)]; % from Kim -- > Pathak
sys.aTu=@(u) [u(2);u(1)];
% ===================  Stationary State  ==================================
sys.X_ss=zeros(6,1);
sys.U_ss=zeros(2,1);
% ===============    Simulation Parameters  ================
% change in SDRE_Simulation.m  for diferent simulations  
sdre.x0=zeros(6,1);     % Initial condition for the system   
sdre.simT=3;            % simulation time
sdre.ts=1e-3;           % simulation step
end
