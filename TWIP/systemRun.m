%% SDRE control design with ROA estimate for a TWIP
% This script runs using 5 additional files:
% systemData.m    :       Solves the optimisation problem and generates ROA 
% systemConfig.m  :       Generates two structures to run simulations
% SDRE_Simulation :       Runs simulations, and saves data for plotting. 
% ---> systemControl.m :  Generates control law inside SDRE_Simulation.m
% ---> dynamics.m      :  Nonlinear system ODE
clc
close all%
clear all
clc
saveImageLatex=0; % Export to tikz
Platform='Windows'; %change to 'Linux' if working on ubuntu, etc.
%% Problem Definition
% set S where the optimisation will be executed
stateLimits=[2, pi/3.0, 4*pi, 5.2*pi];  %systemData
%stateLimits=[2,1, 30*pi/180, pi,179*pi/180, pi]; %systemDataFull
%% Calculation of the Metric
nPoints = [4,4,4,4]; polyGrad = 2; % bad
%nPoints = [2,2,4,4,2,2]; polyGrad = 2; %
BisEnable=1;
Const=0;
[sdre,sys]=systemConfig(stateLimits, nPoints, polyGrad,Const,BisEnable);
%  Load The metric structure

%% Simulation Scenarios
%load('sdreTemp')
clc
%I - Small deviation in Pitch angle with zero velocity vector. 
xI=sys.Tx([0 0 20*pi/180 0 0 0]'); 
%II - Falling(+pitch angle). Forward velocity opposite to driving direction
xII=sys.Tx([0 -.5 35*pi/180 0 15*pi/180 pi/2]'); 
%III - TWIP rapidly falling. Yaw rate opposite to initial yaw deviation
xIII=sys.Tx([0 0 20*pi/180 pi/2 -15*pi/180 -pi]');
%IV - All velocities activated, high yaw rate
xIV=sys.Tx([-1 1 45*pi/180 pi 60*pi/180 3*pi]'); %nice

sdre.fileSim= @WIPfrictionless;
%sdre.fileSim= @WIPfriction;
x0=[xI,xII,xIII,xIV];
sdre.simT=3;
sdre.ts=1e-3;

%sdre.uX= @(x,u) u./50.2810./0.0038.*1.5000; % From torque to voltage
sdre.uX= @(x,u) u.*1; % From torque to voltage
%sdre = rmfield(sdre,'uX');
color=1;
for i=1:3
    sdre.comp=i;
    sdre.x0=x0(:,4);
    if i>2
        color=color+1;
    end
    SDRE_Simulation(sdre,sys,i,1) %1 for WIP System
end

figNames={'TWIPMetric.fig','TWIPNormPosition.fig',...
         'TWIPNormVelocity.fig','TWIPNormInput.fig'};
initValue=2; % in total 6 figures are generated

if saveImageLatex
    for i=1:length(figNames)
        saveas(figure(i+initValue),figNames{i});
    end
    SDRE_SaveTikz('TWIP',0,Platform) %  Export .tikz 
end
%%

figure(3); 
    ylabel('$\lambda(\mathbf{x}),\parallel \mathbf{x} \parallel_M$',...
        'interpreter','latex')
    xlabel('$t$','interpreter','latex')
    %axis([0 finalTime 0 30])
    grid('on')
  
 figure(4); 
    xlabel('$t$','interpreter','latex')
    ylabel('$\parallel \phi \parallel$','interpreter','latex')
    grid('on')
    %axis([0 finalTime 0 30])
  
  figure(5); 
    xlabel('$t$','interpreter','latex')
    ylabel('$\parallel \nu \parallel$','interpreter','latex')
    grid('on')
    %axis([0 finalTime 0 30])
    
  figure(6); 
    xlabel('$t$','interpreter','latex')
    ylabel('$\parallel u \parallel$','interpreter','latex')
    grid('on')
    %axis([0 finalTime 0 30])
