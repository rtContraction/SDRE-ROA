function SDRE_PlotDataSecondOrder(varargin)
%close all
finalPoint=[0;0];
initPoint=[0.8;2.05];
cShade0=[202/255,204/255, 206/255]; % grey
cShade00=[0.6350, 0.0780, 0.1840] 	; 
cShade1=[0, 0.4470, 0.7410];
cShade2=[0.8500, 0.3250, 0.0980];
cShade3=[0.9290, 0.6940, 0.1250];
cShade4=[0.4940, 0.1840, 0.5560]; 
cShade5=[0.4660, 0.6740, 0.1880];

ColorShade={cShade0,cShade1,cShade2,cShade3,cShade4,cShade5};
colorType=varargin{1};
cShade=ColorShade{colorType}; 
start=nargin;
lineWidthColor=1.2;
PointWidthColor=1.5;
% ============= Buffers for the total amount of simulation time============
T_final=[];
x_final=[];
u_final=[];
E_final=[];
dE_final=[];
ur_final=[];
t=0;
out=nargin;
% ============= Buffers for the total amount of simulation time============

for i=start:out
    load(varargin{i},'data'); % reads the data file from each saved file
    T_final=[T_final t+data.tm];
    x_final=[x_final data.x];
    u_final=[u_final data.u];
    ur_final=[ur_final data.ur];
    E_final=[E_final data.E];
    t=t+data.tm(end);           % final array with the time re-arrenged
    if  isfield(data,'dE')
        dE_final=[dE_final data.dE];
    end
end
       % Information of the metric configuration
L=length(T_final);              % Total length of the time vector
n_plots=size(x_final,1);     % dimension of the state vector 
max_fig=4;                      % maximum number of plots per figure
n_fig=round(n_plots/max_fig);   % number of total figures
count=1;

load(varargin{start},'sys'); 
load(varargin{start},'sdre'); 

if sys.Const
    for i=1:n_fig+2
    figure(i);
    hold on
    end
else   
for i=1:n_fig
figure(i);
    for j=1:max_fig
    if count<=n_plots           % check till the maximum number of states
        subplot(max_fig,1,j);
        plot(T_final,x_final(count,:));
        hold on
        %plot(T_final,xr_final(count,:),'m--');
        ylabel(['x',num2str(count)])
        count=count+1;
    end
    end
end
subplot(max_fig,1,3);
    plot(T_final,u_final);
    ylabel('u')
    hold on
    %plot(T_final,ur_final,'m--');
% =============  Contraction Rate ==============================
    
figure(n_fig+1);
    if isequal(sdre.x0,initPoint)
        plot(T_final,ones(length(T_final),1)*sys.cRateChang,'c:',...
            'linewidth',lineWidthColor)
        hold on
        plot(T_final,ones(length(T_final),1)*sdre.cRate,'color',cShade00,...
            'LineStyle','--','linewidth',lineWidthColor)
    end
    %plot(T_final,-dE_final./E_final,'color',cShade);
    if colorType>1
        plot(T_final,-dE_final./E_final,'color',cShade,...
            'linewidth',lineWidthColor);
    else
       plot(T_final,-dE_final./E_final,'color',cShade);
    end
    hold on
% =============  || x ||_P ====================================

 figure(n_fig+2);  
    %plot(T_final,vecnorm(x_final));
    hold on
    a1=sys.lambda_min;
    a2=sys.lambda_max;
    R=sqrt(a2/a1);
    %plot(T_final,sqrt(E_final),'color',cShade);
    
    if colorType>1
        plot(T_final,sqrt(E_final),'color',cShade,...
            'linewidth',lineWidthColor);
    else
       plot(T_final,sqrt(E_final),'color',cShade);
    end
    hold on
    %plot(T_final,vecnorm(x_final).*R.*exp(-sys.cRate.*T_final))

end   
  %Chang Example -Plot Region

figure(n_fig+3);
if isequal(sdre.x0,finalPoint)
    E0=5.45;     
if sys.Const
    xLength=sys.Box(1);
    yLength=sys.Box(2);
    rectangle('Position',[-xLength -yLength 2*xLength 2*yLength],...
        'EdgeColor','c','LineStyle','-.',...
        'linewidth',lineWidthColor)
    hold on
    plot(0,0,':m','linewidth',lineWidthColor)
    plot(-xLength,yLength,'c-.','linewidth',PointWidthColor)
    bracciData([xLength;yLength],lineWidthColor);
    E0=15;
end  
%E0=sdre.E(sdre.x0(1),sdre.x0(2)); %E0=9; %E0=15 for Mconst
[X,Y] = meshgrid(-sys.box(1):0.05:sys.box(1),-sys.box(2):0.05:sys.box(2));
Z=sdre.E(X,Y);

type='k';
if sys.Const
type='k--';
end
contour(X,Y,Z,[E0 E0],type);
else 
    
    if colorType>1
        plot(x_final(1,:),x_final(2,:),'color',cShade,...
            'linewidth',lineWidthColor);
        plot(x_final(1,1),x_final(2,1),'x','color',cShade,...
            'linewidth',PointWidthColor);
    else
        plot(x_final(1,:),x_final(2,:),'color',cShade);
    end
            
end
hold on

end
 