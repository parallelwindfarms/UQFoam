%% Clean-up

clear all;
clc;
close all;

%% Data/ output location

data = getenv('turbChFlowData');
% data = './uqPlots';
caseName = 'n1_FS';
SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/2turbChFlow']);

%% Plot settings

plotp    = 0;
plotU    = 0;
plotUrms = 0;
plotuv   = 0;

plotp    = 1;
plotU    = 1;
plotUrms = 1;
plotuv   = 1;

turbChFlow_main(data, caseName, plotp, plotU, plotUrms, plotuv)