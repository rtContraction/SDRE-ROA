function [Model, Mnum, rDistance,Mdnumf, Coeffs, tbound, feasible, R_nonfeasible,...
          M_nonfeasible, MR_nonfeasible,sys] = systemDataFull(...
          stateLimits, nPoints, polyGrad,Const,BisEnable)
%% 1) create gridded solution of SDRE and split grid into "train"-set 
% that is used to create LMIs for optimization problem and "test"-set to 
% evaluate if conditions for M being a contraction metric also hold  
% outside the chosen training-data
% 1.1) if you want twice as much test-data as training data, use flag 
% "double" ohterwise for test-set to only be slighty shifted in the state
% space region but the same number of points use flag "shifted"
[z_lmi, Acl_pointwise_lmi, Fcl_pointwise_lmi, P_pointwise] = ...
                createLinearCellArrays(stateLimits, nPoints);
%% 2) set up sdpvar-Problem and adapt constraints
% 2.2) define M(z) as symmetric sdpvar with polynomial entries of order 
% sys.met_order, C is the Coefficient-vector of the polynomials, v the 
% vector containing all relevant monomials
n = 6;
z = sdpvar(n,1);  
sys.n=n;
sys.met_order = polyGrad;  
sys.v_ssm = [1 2 3 4 5 6];  
sys.v_ssm_const = []; 
[M,C,v] = SDRE_Metric(z,sys);
tol_clean=1e-12; 
% 2.2) define the closed loop dynmamics
if ~Const
    A_cl = sdpvar(n,n,'full');  
    f_cl = A_cl*z; 
% 2.3) define M_dot: time derivative of contraction metric M(x)
    Md = SDRE_Md(n, z, f_cl, M); 
else
    Md=0*eye(n,n);
end
% 2.4) define t as coverge rate
% t can be chosen constant for faster solving or maximized using bisection
if BisEnable
t = sdpvar(1);
else
t=1;
end
% 2.5) add contraints  
% constraints can be for M-P being positive definite or relaxed to M-eps*I
% being positive definte, where P is the local ARE-Solution and I is the 
% identiy matrix
Model = setConstraints(z_lmi, Acl_pointwise_lmi, Fcl_pointwise_lmi,...
                           P_pointwise, M, Md, A_cl, z, t,Const);   
                       
%% 3) start optimization - used solver is mosek
% 3.1)  dependent on how t was chosen different optimier-calls need to be
% made
ops = sdpsettings('solver', 'mosek');
Objective=[];
if BisEnable
    ops.bisection.absgaptol=1e-4; % tolerance to reduce the gap
    Objective = -t;
    diagnostics = bisection(Model, Objective, ops);
    tbound = value(t)*0.9;         
    if tbound < 0.0
        error('Negative contraction rate, simulation aborted');
        return
    end
else
    tbound=t;
    diagnostics = optimize(Model,Objective,ops);
end

% 3.2) postprocess result from optimization 

Coeffs = value(C);
Mnum = replace(M, C, Coeffs);
Mnum=clean(Mnum,tol_clean);  

rDistance=matlabFunction(SDRE_SDP2Sym(z,z.'*Mnum*z),'Vars',...
{SDRE_SDP2Sym(z,z(1)),SDRE_SDP2Sym(z,z(2)),...
SDRE_SDP2Sym(z,z(3)),SDRE_SDP2Sym(z,z(4)),...
SDRE_SDP2Sym(z,z(5)),SDRE_SDP2Sym(z,z(6))}); % norm(x)_M^2 = x^T M(x) * x
                                             % use to plot ROA
[~,~,Mnum] = SDRE_Posprocessing...
    (Mnum, v, z, tol_clean, sys.v_ssm);     % Function handle M(x)

if ~Const
    Mdnum=replace(Md,C,Coeffs);             % Solution dM(x)
    Mdnum=clean(Mdnum,tol_clean);
    Mdnumf=matlabFunction(SDRE_SDP2Sym([z;A_cl(:)],Mdnum),'Vars',...
    {SDRE_SDP2Sym([z;A_cl(:)],[z;A_cl(:)])}); % Function handle dM(x)
else
    Mdnumf=@(x) zeros(n,n);
end

%% 4) validate metric M that was determined through optimizer
% 4.1) check if conditions are actually fulfilled in train-set - can be
% omitted as it has never failed if solver converged
% [feasible, M_nonfeasible, R_nonfeasible, MR_nonfeasible]= ...
%  evaluateMetrics(z_lmi, Acl_pointwise_lmi, Fcl_pointwise_lmi, Mnum, Md, A_cl, z, tbound);
feasible=1;M_nonfeasible=1;R_nonfeasible=1; MR_nonfeasible=1;
end

function [z_lmi, Acl_pointwise_lmi, Fcl_pointwise_lmi, P_pointwise_lmi] =...
          createLinearCellArrays(stateRange, nPoints)
      
    %% 1) first, assure that the number of points is uneven for all states 
    %% such that x = 0 is always included in evaluation 
    nPoints(mod(nPoints,2) == 0) = nPoints(mod(nPoints,2) == 0) + 1;
     if ~isempty(nPoints(mod(nPoints,2) == 0))
         fprintf("Even number in number of points!\n")
     end
    
    %% 2) first create sets for LMIs in optimization problem
    %% 2.1) define range and initialize empty arrays
    x_steps_lmi = nPoints(1);
    xdot_steps_lmi = nPoints(2);
    theta_steps_lmi = nPoints(3);
    thetad_steps_lmi = nPoints(4);
    psi_steps_lmi = nPoints(5);
    psid_steps_lmi = nPoints(6);
    
    x_lmi = linspace(-stateRange(1), stateRange(1), nPoints(1));
    xdot_lmi = linspace(-stateRange(2), stateRange(2), nPoints(2));
    theta_lmi = linspace(-stateRange(3), stateRange(3), nPoints(3));
    thetadot_lmi = linspace(-stateRange(4), stateRange(4), nPoints(4));
    psi_lmi = linspace(-stateRange(5), stateRange(5), nPoints(5));
    psidot_lmi = linspace(-stateRange(6), stateRange(6), nPoints(6));   
    
    p_total_lmi = nPoints(1)*nPoints(2)*nPoints(3)*nPoints(4)...
        *nPoints(5)*nPoints(6);
    
    z_lmi = cell(nPoints(1), nPoints(2), nPoints(3), nPoints(4), ...
        nPoints(5), nPoints(6));
    Acl_pointwise_lmi = cell(nPoints(1), nPoints(2), nPoints(3), nPoints(4), ...
        nPoints(5), nPoints(6));
    Fcl_pointwise_lmi = cell(nPoints(1), nPoints(2), nPoints(3), nPoints(4), ...
        nPoints(5), nPoints(6));      
    P_pointwise_lmi = cell(nPoints(1), nPoints(2), nPoints(3), nPoints(4), ...
        nPoints(5), nPoints(6));    
     %% 2.2) Pointwisely solve SDRE and compute required matrices
    barGUI = waitbar(1/p_total_lmi, 'Initializing...', 'Name', ...
                    'SDRE Simulation running - Train set');
    tic
    
    for h = 1:x_steps_lmi
        for i = 1:xdot_steps_lmi
            for j = 1:theta_steps_lmi
                for k = 1:thetad_steps_lmi
                    for l = 1:psi_steps_lmi
                        for m = 1:psid_steps_lmi
                       

                %% SDRE - pointwise evaluation and computation of generalised
                %% gradient F
                z_curr = [x_lmi(h), xdot_lmi(i), theta_lmi(j), ...
                          thetadot_lmi(k), psi_lmi(l), psidot_lmi(m)]';
                [A, B] = systemSDC(z_curr);
                [Q, R] = systemWeightingMatrices(z_curr);   
                dA_dx_x = systemJacobianOpenLoop(z_curr);
                [P, K, ~, ~] = icare(A,B,Q,R,[],[],[]);
                dBKx_dx = systemJacobianStateFeedback(theta_lmi(j), ...
                                                 K(1,3), K(2,3));                            
                Acl_pointwise_lmi{h,i,j,k,l,m} = A - B * K;
                Fcl_pointwise_lmi{h,i,j,k,l,m} = A + dA_dx_x - B*K - dBKx_dx;
                P_pointwise_lmi{h,i,j,k,l,m} = P;
                z_lmi{h,i,j,k,l,m} = z_curr;
                
                curr_point = m+(l-1)*psid_steps_lmi+...
                    (k-1)*psid_steps_lmi*psi_steps_lmi+...
                    (j-1)*psid_steps_lmi*psi_steps_lmi*thetad_steps_lmi+...
                    (i-1)*psid_steps_lmi*psi_steps_lmi*thetad_steps_lmi*theta_steps_lmi+...
                    (h-1)*psid_steps_lmi*psi_steps_lmi*thetad_steps_lmi*...
                    theta_steps_lmi*xdot_steps_lmi;        
                waitbar(curr_point/p_total_lmi, barGUI, sprintf('Simulation running'));
                        end
                    end
                end
            end
        end
    end
    toc
    waitbar(1, barGUI, 'Finishing');
    close(barGUI);           
end

function Model = setConstraints(z_Array, Acl_Array, Fcl_Array, P, M,...
                                     Md, A_cl, z, t,Const)

%% create empty Model to append pointwise LMIs
Model = []; 
n = 6;
x_steps = size(Fcl_Array,1);
v_steps = size(Fcl_Array,2);
theta_steps = size(Fcl_Array,3);
thetad_steps = size(Fcl_Array,4);
psi_steps = size(Fcl_Array,5);
psid_steps = size(Fcl_Array,6);

dM=jacobian(M,z);
%% create waitbar to update about solver progress
p_lmi = x_steps*v_steps*theta_steps*thetad_steps*psid_steps*psi_steps; %total number of points for LMIs
barGUI = waitbar(1/p_lmi, 'Initializing...', 'Name', 'Constraints');

%% set LMIs either for M to be positive definite or M-P being positive 
%% definite. Either way also append condition for M to be a contration
%% metric
tic
    for h = 1:x_steps
        for i = 1:v_steps
            for j = 1:theta_steps
                for k = 1:thetad_steps
                    for l = 1:psi_steps
                        for m = 1:psid_steps
           z_curr = z_Array{h,i,j,k,l,m};
           F_curr = Fcl_Array{h,i,j,k,l,m};
           Acl_curr = Acl_Array{h,i,j,k,l,m};
            %% add postive semi-definite P as pointwise LMI for M(x) to Model
           P_curr = P{h,i,j,k,l,m};
           Model = [Model,1e5*eye(n,n)>=replace(M,z,z_curr)>= 1e-5*eye(n,n) ]; %e4,e-5 
           if isequal(z_curr, zeros(n,1))
               z_curr
               %Model = [Model, norm(replace(M,z,z_curr)-P_curr,2)<=tol_norm ];
               %Model = [Model, replace(M, z, z_curr) == P_curr]
%                states=[3];
%                for s=1:length(states)     
%                Model = [Model, replace(dM(:,:,states(s)),z,z_curr) == zeros(n,n) ];
%                end
           end
               % else Model = [Model, replace(M, z, z_curr) >= P_curr];     
           %% optional: set upper bound for M(x)
           % Model = [Model, replace(M, z, z_curr) <= 1e7.*eye(6)]; 
           
           %% append condition for M(x) to be contraction metric
           % ===========  Inequality condition  for Contraction  =========
           % ===========  M_dot + M*F + F'*M < -beta *M  =================
            if ~Const
                 Mdi = replace(Md, A_cl(:), Acl_curr(:)); 
            else 
                Mdi=0*eye(n,n);
            end  
            Model = [Model, replace((transpose(F_curr)*M + M*F_curr + ...
                       Mdi + (2*t)*M), z, z_curr) <= 0];
            %% update waitbar           
               curr_point = m+(l-1)*psid_steps+...
                    (k-1)*psid_steps*psi_steps+...
                    (j-1)*psid_steps*psi_steps*thetad_steps+...
                    (i-1)*psid_steps*psi_steps*thetad_steps*theta_steps+...
                    (h-1)*psid_steps*psi_steps*thetad_steps*...
                    theta_steps*v_steps;              
            waitbar(curr_point/p_lmi, barGUI, sprintf('Adapting Contraints'));
                        end
                    end
                end
            end
        end
    end
toc
waitbar(1, barGUI, 'Finishing');
close(barGUI);

end

function [feasible, M_nonfeasible, R_nonfeasible, MR_nonfeasible] = ...
          evaluateMetrics(z_Array, Acl_Array, Fcl_Array, Mnum, Md, A_cl,...
                          z, t)
%% initialize empty arrays
feasible_counter = 0;
M_nonfeasible_counter = 0;
MR_nonfeasible_counter = 0;
R_nonfeasible_counter = 0;
feasible = [];
M_nonfeasible = [];
R_nonfeasible = [];
MR_nonfeasible = [];

%% detmermine array sizes
xdot_steps = size(Fcl_Array,1);
theta_steps = size(Fcl_Array,2);
thetad_steps = size(Fcl_Array,3);
psid_steps = size(Fcl_Array,4);

%% create waitbar to update on progress
p_total = xdot_steps*theta_steps*thetad_steps*psid_steps;
barGUI = waitbar(1/p_total, 'Initializing...', 'Name', 'Evaluation of Metric');
    
tic
for i = 1:xdot_steps
    for j = 1:theta_steps
        for k = 1:thetad_steps
           for l = 1:psid_steps
          
           %% "double"-test-set might be empty for some entries
           if isempty(z_Array{i,j,k,l})
               continue
           end
           
           %% evaluate M(z) obtained from optimization and see if 
           %% conditions on contraction are given
           z_curr = z_Array{i,j,k,l};
           Acli = Acl_Array{i,j,k,l}; % 
           Mdi = replace(Md, A_cl(:), Acli(:));   %
           Mdi = replace(Mdi, z, z_curr); 
           R = double(Mdi) + Mnum(z_curr)*Fcl_Array{i,j,k,l} + ...
               Fcl_Array{i,j,k,l}'*Mnum(z_curr) + 2*t*Mnum(z_curr);
           
           %% update waitbar 
           curr_point = l+(k-1)*psid_steps+...
                    (j-1)*psid_steps*thetad_steps+...
                    (i-1)*psid_steps*thetad_steps*theta_steps  ; 
           waitbar(curr_point/p_total, barGUI, ...
                   sprintf('min/max Eigenvalues of M(x)/R(x): %0.1f, %0.1f',...
                   min(real(eig(Mnum(z_curr)))), max(real(eig(R))) ));

            if max(real(eig(R))) <=0 && min(real(eig(Mnum(z_curr)))) > 0
%                z_curr
                feasible_counter = feasible_counter + 1;
                feasible = [feasible, z_curr];
            elseif max(real(eig(R)))>0 && min(real(eig(Mnum(z_curr)))) > 0
%                fprintf('-R(x) is not positive semi-definite at:\n')
                %z_curr
                R_nonfeasible_counter = R_nonfeasible_counter + 1;
                R_nonfeasible = [R_nonfeasible, z_curr];
            elseif max(real(eig(R)))<=0 && min(real(eig(Mnum(z_curr))))<=0
%                fprintf('M(x) is not postive definite at:\n')
%                    z_curr
                M_nonfeasible_counter = M_nonfeasible_counter + 1;
                M_nonfeasible = [M_nonfeasible, z_curr];
            else
%                fprintf('Neither M(x) nor -R(x) positive definite at:\n')
%                     z_curr;
                MR_nonfeasible_counter = MR_nonfeasible_counter + 1;
                MR_nonfeasible = [MR_nonfeasible, z_curr];
            end
                
           end
        end
    end
end

waitbar(1, barGUI, 'Finishing');
close(barGUI);

%% summarize how many points are not feasible/ feasible
fprintf('Of a total of %d checked points in %d conditions for M(x) to be\n',...
        p_total, feasible_counter) 
fprintf('a contraction metric are fulfilled, and in %d not.\n',...
        M_nonfeasible_counter + R_nonfeasible_counter + MR_nonfeasible_counter)
pause(5)
end





