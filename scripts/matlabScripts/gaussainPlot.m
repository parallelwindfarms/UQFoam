%% Clean
 
clear all;
clc;
close all;

%% Gaussian

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/matlabScripts/export_fig/export_fig-master']);

x = linspace(-3,3,3000);
y = normpdf(x,0,1);
purple = [69,2,86]/255;
yellow = [249,231,33]/255;
green = [175, 231,33]/255;
patch(x,y, green,'FaceAlpha', 0.75,'EdgeColor','k');
set(gca,'visible','off')
set(gca,'xtick',[])
set(gca,'XColor', 'none','YColor','none')

% PNG
% location = '/home/jigar/Dropbox/1Groningen/PhD/14PaperThesis/1_UQ_IPC_SIAM/figures';
% figname = [location '/gaussian_yellow'];
export_fig gaussian_yellow -transparent -r1200 -png
% export_fig(figname,'-transparent','-r1200','-png'); 
