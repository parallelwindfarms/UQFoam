%% Clean-up

clear all;
clc;
close all;

%% Data/ output location

data = getenv('hpFlowData');
caseName = '1hpFlow';

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/1hpFlow']);

%% Plot settings

plotp0_c = 0;
plotU0_c = 0;
plotGradp = 0;
plotU_xByd = 0;
plotU_xByd_offset = 0;

plotp0_c = 1;
plotU0_c = 1;
plotGradp = 1;
plotU_xByd = 1;
plotU_xByd_offset = 1;

hpFlow_main(data, caseName, plotp0_c, plotU0_c, plotGradp, plotU_xByd, plotU_xByd_offset)

