%% Clean-up

clear all;
clc;
close all;
restoredefaultpath;

%% Data/ output location

data = getenv('periodicHillData');
% data = './uqPlots';
caseName = '2D_0DET_1baseline_kEpsilon';

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/7periodicHill']);

%% Plot settings

plotU_xByd    = 0;
plotUU_xByd   = 0;
plotVV_xByd   = 0;
plotUV_xByd   = 0;

plotU_xByd    = 1;
plotUU_xByd   = 1;
plotVV_xByd   = 1;
plotUV_xByd   = 1;

periodicHill_main_det(data, caseName, plotU_xByd, plotUU_xByd, plotVV_xByd, plotUV_xByd)