function [Mnumf,rDistance,Mdnumf,trajBound,sys] = systemData...
                                              (xrange,nPoints,degree,Const)
% Second Order System as in Bracci-2006 and Chang-2019
%==================  Optimasation Tolerances  =============================
tol_clean=1e-6; %tolantence e-12
tol_norm=1e-2;
%==================   System Properties  ==================================
b1=0;
b2=1;
B=[b1;b2];
n=2;
Jf = @(x) [2*x(1)  1;0 0];
A_x = @(x) [x(1) 1; 0 0];
%==================   SDRE Parameters =====================================
R=1;
Q = diag([100,100]);
%==================  Gridding Terms =======================================
xPoints1=nPoints(1)+1; 
xPoints2=nPoints(2)+1;
xP1=linspace(-xrange(1),xrange(1),xPoints1);
xP2=linspace(-xrange(2),xrange(2),xPoints2);
Acl=cell(xPoints1,xPoints2);
Fcl=cell(xPoints1,xPoints2);
L=xPoints1*xPoints2; % Total number of points
% ===========================  SDP Variables  =============================
z=sdpvar(n,1); % free variable for optimasation
beta=sdpvar(1);
sys.n=n;
sys.met_order = degree;  
sys.v_ssm = [1 2];   % selection of the states x1,x2 
sys.v_ssm_const=[];  % partitions which might be constant
[M,C,v]=SDRE_Metric(z,sys); % metric, decision variables, Monomials  
if ~Const
    c_l=sdpvar(n*n,1,'full'); 
    f_cl=[c_l(1) c_l(3);c_l(2) c_l(4)]*z; % f_cl = A_cl(x)*x
    Md=SDRE_Md(n,z,f_cl,M);                % Mdot = dM/dx * f_cl
else
    Md=0*eye(n,n);                        % assume Mdot= 0
end
dM=jacobian(M,z);                         % tensor dM(x)/dx_i, i=1...n 
% ===================  1. Gridding Points  ================================
Model=[];
barGUI = waitbar(1/L,'1','Name','SDRE Simulation running');               
for j=1:xPoints2
    for i=1:xPoints1
        x=[xP1(i);xP2(j)];
        %  ==================== SDRE ======================================
        [Psol,gainSol,eigSol]=icare(A_x(x),B,Q,R);  %Q(x3) ?, R(x3)?
        Acli=A_x(x) - B*gainSol;
        Acl{i,j}=Acli;
        %  ==================== YALMIP SDP-Model ==========================
        if isequal(x,zeros(n,1)) % flatting dM/dx at x=0
            Model = [Model, replace(dM(:,:,1),z,x) == zeros(n,n) ];
            Model = [Model, replace(dM(:,:,2),z,x) == zeros(n,n) ];                   
        end
        Model = [Model,1.5e2*eye(n,n) >=replace(M,z,x)>= 1e-1*eye(n,n) ]; 
        % Model = [Model,1.4*Psol >= replace(M,z,x)>= .6*Psol ]; % M = P  
        % ===========  Contraction Condition  =============================
        if ~Const
            Mdi=replace(Md,c_l,Acli(:)); % replace current x in dM(x) 
        else 
            Mdi=0*eye(n,n);
        end    
        F=Jf(x)-B*gainSol; % closed-loop system
        Fcl{i,j}=F;
        H= M*F + F.'*M;
        Rcl = -(Mdi + H);
        Model = [Model,replace(Rcl - beta*M,z,x)>=0];
        waitbar((i+xPoints2*(j-1))/L,barGUI,...
            sprintf('Pointwise evaluation at:%0.5f,%0.5f',x(1),x(2)))
    end
end
waitbar(1,barGUI,'Solving Optimisation...');
Objective = -beta;
ops = sdpsettings('solver', 'mosek');
if ~Const
    ops.bisection.absgaptol= 1e-4; % gap tolerance for  t
else
    ops.bisection.absgaptol= 1e-5; % gap tolerance for  t
end
diagnostics = bisection(Model,Objective,ops); 
% ===========================  Posprocessing  =============================                 
Mnum=replace(M,C,value(C)); % solution of M(x)
Mnum=clean(Mnum,tol_clean);                 
rDistance=matlabFunction(SDRE_SDP2Sym(z,z.'*Mnum*z),'Vars',...
{SDRE_SDP2Sym(z,z(1)),SDRE_SDP2Sym(z,z(2))});  % norm(x)_M^2 = x^T M(x) * x
                                             % use to plot ROA
if ~Const
    Mdnum=replace(Md,C,value(C));            % Solution dM(x)
    Mdnum=clean(Mdnum,tol_clean);
    Mdnumf=matlabFunction(SDRE_SDP2Sym([z;c_l],Mdnum),'Vars',...
    {SDRE_SDP2Sym([z;c_l],[z;c_l])});         % Function handle dM(x)
else
    Mdnumf=@(x) [0,0;0,0];
end
[~,~,Mnumf]= SDRE_Posprocessing...
            (Mnum,v,z,tol_clean,sys.v_ssm); % Function handle M(x)
% ================  Second Part, check that M(x) > 0 and -R(x)>0  =========
betaBound = value(beta)*0.9;
M0=Mnumf(zeros(2,1));
lambda_min=min(real(eig(M0))); % Capture  a initial stimate of the 
lambda_max=max(real(eig(M0))); % slowest and fastest eigenvalue of M(x)
s=0;                           % Counting element                            
for j=1:xPoints2
    for i=1:xPoints1
        x=[xP1(i);xP2(j)];
        Acli=Acl{i,j};
        if ~Const
            Mdi=replace(Mdnum,c_l,Acli(:));  
        else
            Mdi= 0*eye(n,n);
        end
        H= Mnum*Fcl{i,j} + Fcl{i,j}.'*Mnum;
        R=Mdi+ H + 2*betaBound*Mnum;
        Rx=double(replace(R,z,x));
        Mx=double(replace(Mnum,z,x));
        if min(real(eig(-Rx)))>=0 && min(real(eig(Mx)))>0
            if min(real(eig(Mx)))<lambda_min
               lambda_min=min(real(eig(Mx)));   % new slowest eig
            end   
            if max(real(eig(Mx)))>lambda_max
               lambda_max=max(real(eig(Mx))); % new fastest eig
            end
        else
            s=s+1;% Points where the condition was not satisfied
        end
        
    end
end
trajBound=[betaBound;lambda_min;lambda_max];
waitbar(1,barGUI,sprintf('No solution in %0.5f points',s));
waitbar(1,barGUI,'Closing');
close(barGUI)
end
