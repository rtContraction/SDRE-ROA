%% SDRE control design with ROA estimate for a second order system
% This script runs using 5 additional files:

% systemData.m    :       Solve optimisation problem and compute ROA 
% systemConfig.m  :       Generate two structures to run simulations
% SDRE_Simulation :       Run simulations  and save data for plotting. 
% ---> systemControl.m :  Generate control law inside SDRE_Simulation.m
% ---> dynamics.m      :  Nonlinear system ODE
close all%
clear all
clc
saveImageLatex=0; % Export to tikz
Platform='Windows'; %change to 'Linux' if working on ubuntu, etc.
%% Problem Definition
xPoints=20;yPoints=20; % set S where the optimisation will be executed
%% Calculation of the Metric
% sdre structure contain information of the metric and controler
% sys contain information of the nonlinear system
% both are used to run  Simulations  and data processing for plotting

% Case 1, Nonlinear Metric with Mdot(x) = 0
xLength=1;yLength=1;% This is the maximum range where a solution is found
metricDegree=4;Mconst=1;
[sdre0,sys0]=systemConfig([xLength;yLength],[xPoints;yPoints],...
    metricDegree,Mconst);

% Case 2, Nonlinear Metric with Mdot(x)
xPoints=30;yPoints=30;
xLength=5;yLength=5;% 8 Metric for control
metricDegree=2; Mconst=0;
[sdre1,sys1]=systemConfig([xLength;yLength],[xPoints;yPoints],...
 metricDegree,Mconst);
%% Selection of Trajectories 
% Collection of trajectories to be plotted
X0=[0;0];
% A(I)-Quadrant
A1=[0.8;2.05];A2=[.25;3.15];A3=[0;3.44];
xA=[A1,A2,A3];
% B(II)-Quadrant
B1=[-1.41;3.95];B2=[-1.25;1.95];
xB=[B1,B2];
% C(III)-Quadrant
C1=[-0.75;-2.2];C2=[-0.15;-3.25];
xC=[C1,C2];
% D(IV)-Quadrant
D1=[0.75;-4.35];D2=[1.35;-3.1]; 
xD=[D1,D2];
% Outside
O1=[0;-5]; % to colorized
O2=[-5;0];O3=[-2;-5];
xO=[O1,O2,O3];
% Colored 
P1=[-1.3;4.95];     % Blue
P2=[1.3;-4.95];   % Orange
P3=[1.2;1.2];       % Yellow
P4=[-5;5];          % Purple
P5=[2;2];           % Green
xP=[P1,P2,P3,P4,P5];

x0=[xA,xB,xC,xD,xO,X0,xP];  % Collection of all trajectories
intLength=size(x0,2)-size(xP,2); % Trajectories display in grey
%% Simulation 
sdre1.simT=5;  % Simulation time
sdre1.ts=1e-2; % Sample Time
sdre1.comp=2;  % 1= CLF using Backstepping, 2 = SDRE + Contraction Theory
color=1;       % Flag to plot using color
sys0.Box=sys1.box;
SDRE_Simulation(sdre0,sys0,color,0) % 0 = This Example , 1 = TWIP
for i=1:size(x0,2)
    sdre1.x0=x0(:,i); % Change the initial condition
    if i>intLength
        color=color+1;    % Plot with color line
    end
    SDRE_Simulation(sdre1,sys1,color,0) % Last run to plot ROA
end

% JUST FOR US TO CREATE TIKZ. no for actual repo (line 73 -81)
figNames={'states.fig','rate.fig','norm.fig','region.fig'};
for i=1:length(figNames)
saveas(figure(i),figNames{i});
end

if saveImageLatex
    SDRE_SaveTikz('secondOrder',intLength,Platform) %  Export .tikz 
end
%% Formatting the plots  -- TO DO (only for MATLAB display)

finalTime=5;
figure(1);
figure(2);
ylabel('$\lambda(\mathbf{x})$','interpreter','latex')
xlabel('$t$','interpreter','latex')
axis([0 finalTime 0 30])
grid('on')
    
figure(3); 
xlabel('$t$','interpreter','latex')
ylabel('$\parallel \mathbf{x} \parallel_M$','interpreter','latex')
grid('on')
axis([0 finalTime 0 30])
    
figure(4);
axis([-xLength xLength -yLength yLength])
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')




% 
