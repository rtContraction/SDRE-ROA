function SDRE_PlotDataWIP(varargin)
%close all
cShade0=[128/255,128/255, 128/255]; % grey
cShade1=[0.6350, 0.0780, 0.1840] 	; 
cShade2=[0, 0.4470, 0.7410];
cShade3=[0.8500, 0.3250, 0.0980]; 
cShade4=[0.4940, 0.1840, 0.5560]; 
cShade5=[0.9290, 0.6940, 0.1250];
ColorShadeA={cShade2,cShade3,cShade5};
ColorShadeB={cShade2,cShade3,cShade5};
cShadeA=ColorShadeA{varargin{1}};
cShadeB=ColorShadeB{varargin{1}}; 
start=nargin;
lineWidthColor=1.2;
% ============= Buffers for the total amount of simulation time============
T_final=[];
x_final=[];
u_final=[];
E_final=[];
dE_final=[];
xr_final=[];
ur_final=[];
t=0;
out=nargin;
% ============= Buffers for the total amount of simulation time============

for i=start:out
    load(varargin{i},'data'); % reads the data file from each saved file
    T_final=[T_final t+data.tm];
    x_final=[x_final data.x];
    u_final=[u_final data.u];
    xr_final=[xr_final data.xr];
    ur_final=[ur_final data.ur];
    E_final=[E_final data.E];
    t=t+data.tm(end);           % final array with the time re-arrenged
    if  isfield(data,'dE')
        dE_final=[dE_final data.dE];
    end
end
load(varargin{start},'sys'); 
load(varargin{start},'sdre');        % Information of the metric configuration
L=length(T_final);              % Total length of the time vector
n_plots=size(x_final,1);     % dimension of the state vector 
max_fig=4;                      % maximum number of plots per figure
n_fig=round(n_plots/max_fig);   % number of total figures

count=1;
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
    if sdre.comp==3
        plot(T_final,ones(length(T_final),1)*2*sys.cRate,'color',cShade2,...
            'LineStyle','--','linewidth',lineWidthColor);
        hold on
        plot(T_final,-dE_final./E_final,'color',cShade3,...
            'linewidth',lineWidthColor);
%        plot(T_final(1:end-1),diff(log(E_final)),'color',cShade3,...
%             'linewidth',lineWidthColor);
%         plot(T_final(1:end-1),diff(E_final)./E_final(1:end-1),'color',cShade3,...
%             'linewidth',lineWidthColor);
        plot(T_final,sqrt(E_final),'color',cShade5,...
            'linewidth',lineWidthColor);
    end
% =============  x,v plots ====================================
    figure(n_fig+2);
    x_norm=[x_final(1,:);x_final(3,:);x_final(5,:)];
    plot(T_final,vecnorm(x_norm),'color',cShadeA,'linewidth',lineWidthColor);
    hold on
    figure(n_fig+3);
    v_norm=[x_final(2,:);x_final(4,:);x_final(6,:)];
    plot(T_final,vecnorm(v_norm),'color',cShadeA,'linewidth',lineWidthColor);
    hold on    
% =============  inputs ====================================
%     figure(n_fig+3);
%     plot(T_final,u_final(1,:),'color',cShadeA);
%     hold on
%     plot(T_final,u_final(2,:),'color',cShadeB);
    figure(n_fig+4);
    plot(T_final,vecnorm(u_final),'color',cShadeA,'linewidth',lineWidthColor);
    hold on
    
end
 