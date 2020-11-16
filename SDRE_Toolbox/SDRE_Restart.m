function SDRE_Restart(varargin)
% function to reset all values of the ccm_config files before executing
if nargin > 0
    set(0,'DefaultFigureWindowStyle','normal'); % one window per figure
    clear all                                   % clean also the cache
    clc 
    close all
else
    clear variables                     % clean workspace variables
    set(0,'DefaultFigureWindowStyle','docked'); 
    yalmip('clear')                     % clean LMI variables
    yalmip('clearsos')                  % clean SoS variables
    yalmip('clearsolution') %clear the variables in F                                    
    disp('clean completed') % clean main display
end
end
% checked by j.perez 29.10.2020
