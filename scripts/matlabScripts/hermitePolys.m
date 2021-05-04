%% Clean
 
clear all;
clc;
close all;

%% Hermite Polynomials

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/matlabScripts/export_fig/export_fig-master']);

fSize = 30;
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', fSize);
% set(0,'DefaultLegendFontSize',fSize);
% set(0,'DefaultTextFontSize',fSize);

% set(groot, 'units', 'inches', 'position', [0 0 8 4])
set(groot, 'defaultFigureUnits','inches')
set(groot, 'defaultFigurePosition',[2.5 1.5 8 6])

% set(groot, 'defaultFigurePaperPositionMode', 'manual');
set(groot, 'defaultFigurePaperUnits', 'inches');
set(groot, 'defaultFigurePaperPosition', [2.5 1.5 8 6]);

syms x y
fplot(hermiteH(0:4,x),'LineWidth',2)
axis([-2 2 -30 30])
yticks([-30 0 30]);
grid off

leg1 = legend('$\phi_0(x)$', '$\phi_1(x)$', '$\phi_2(x)$', '$\phi_3(x)$', 'Location', 'south');
set(leg1,'Visible','on');
set(leg1,'Color','none');
set(leg1,'EdgeColor','none');

% PNG
set(gca, 'Color', 'none');
location = '/home/jigar/Dropbox/1Groningen/PhD/14PaperThesis/1_UQ_IPC_SIAM/figures';
figname = [location '/hermitePolynomials'];
% export_fig hermitePolynomials -transparent -r1200 -png
export_fig(figname,'-transparent','-r1200','-png'); 