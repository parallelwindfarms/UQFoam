%% Clean-up

clear all;
clc;
close all;

%% Data/ output location

data = getenv('flowOvrCylData');
% data = './uqPlots';
caseName = 'detCase_M1';

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/3flowOvrCyl']);

%% Plot settings

plotU_c       = 0;
plotU_xByd    = 0;
plotV_xByd    = 0;
plotUU_xByd   = 0;
plotVV_xByd   = 0;
plotCd        = 0;
plotCl        = 0;
plotCp        = 0;

plotU_c       = 1;
plotU_xByd    = 1;
plotV_xByd    = 1;
plotUU_xByd   = 1;
plotVV_xByd   = 1;
% plotCd        = 1;
% plotCl        = 1;
% plotCp        = 1;

flowOvrCyl_main_det(data, caseName, plotU_c, plotU_xByd, plotV_xByd,...
                    plotU_xByd, plotVV_xByd, plotCd, plotCl, plotCp)