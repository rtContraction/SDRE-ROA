% =======================  Control synthesis   ============================
function SDRE_Simulation(varargin)
% =================  Path to the configuration file =======================
sdre=varargin{1};
sys=varargin{2};
color=varargin{3};
plant=varargin{4};
n=sdre.n;
m=sdre.m;
% ===========================  Data Settings  =============================
T=sdre.simT;  
x=sdre.x0; % stablish the initial condition
% =======================  Data Storage  ==================================
ts=sdre.ts; tm=0:ts:T; L=length(tm); data.tm=tm;   % time vector - t
if  isfield(sdre,'set')
    data.xr=sys.X_ss; 
else
    data.xr=repmat(sys.X_ss,1,L);
end
data.ur=repmat(sys.U_ss,1,L);                     % pair (xr,ur) -Desired
data.x=zeros(n,L); data.u=zeros(m,L);             % pair (x,u)   -Plant
data.E=zeros(1,L);          
data.dE=zeros(1,L);
% =======================  Start from a previous simulation ===============
xr = data.xr(:,1);                           % takes xr as x_ss                      
ur = data.ur(:,1);                          % takes ur as u_ss
%% Simulation
for t=1:L
    % ===================  Feedforwd tracking  ============================
    if sdre.track                                                           
        ur = sdre.ur(data.tm(t),xr);
        data.ur(:,t) = ur;
        data.xr(:,t) = xr;
    end
    % =======================  Transformation  ============================
    xr=sys.Tx(data.xr(:,t));
    x=sys.Tx(x);
    ur=sys.Tu(ur);
    % ====================  Control  ======================================
    [u,E,dE]=systemControl(x,xr,ur,sdre,sys);
    % =================== Antitransformations  ============================
    x=sys.aTx(x);
    u=sys.aTu(u);
    xr=sys.aTx(xr);
    % ================= Save Data Antitransformations  ====================
    data.x(:,t)=x;
    data.u(:,t)=u;
    data.E(t)=E;
    data.dE(t)=dE;
    if isfield(sdre,'uX')
        data.u(:,t)=sdre.uX(x,sys.Tu(u));
        data.x(:,t)=sys.Tx(x);
    end
    if abs(E) > 1e20
        f = waitbar(1/L,'1','Name','SDRE Simulation running');
        waitbar(L/L,f,...
            sprintf('Unbounded Energy. %0.2f',E))
            waitbar(1,f,'Finishing');
            disp("======= ERROR at Time =======");
            disp([t,tm(t)]); % Time when the error ocurres 
            disp(E); % Time when the error ocurres 
            data.E(t:end)=E; % for polting a maximum value
            data.x(:,t:end)=repmat(data.x(:,t-1),1,L-t+1); %last state
        waitbar(1,f,'Finishing');
        close(f)
        break
    else
    % =========== ODE Solver for the next step ============================
        if  isfield(sdre,'fileSim')
            x = SDRE_NextState(ts,sdre.fileSim,x,u);
            if sdre.track
                xr = SDRE_NextState(ts,sdre.fileSim,xr,ur);
            end
        else
            x = SDRE_NextState(ts,@dynamics,x,u);
            if sdre.track
                xr = SDRE_NextState(ts,@dynamics,xr,ur);
            end
        end
        
    end
end
% ======================== Save Temporary Files ===========================
save('sdreTemp.mat','sys','sdre','data')
if plant==0 %% Chang2009 Example
    SDRE_PlotDataSecondOrder(color,'sdreTemp.mat')
else
    SDRE_PlotDataWIP(color,'sdreTemp.mat')
end
end


