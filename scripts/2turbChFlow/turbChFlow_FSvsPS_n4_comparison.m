%% Clean-up

clear all;
clc;
close all;

data = getenv('turbChFlowData');

%% Load data and scripts

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/matlabScripts']);
addpath([SCRIPTS '/matlabScripts/matlab2tikz-master/src']);
addpath([SCRIPTS '/2turbChFlow']);

load([data '/MC_data.mat']);

order       = 4;
numCases    = 3;
UMean       = zeros(52,numCases);
UMeanSigma  = zeros(52,numCases);
RMean       = zeros(52,numCases,4);
RMeanSigma  = zeros(52,numCases,4);

for n = 1:numCases

if n == 1
    load([data '/n'  num2str(order) '_FS/n' num2str(order) '_FS.mat']);
elseif n == 2
    load([data '/n'  num2str(order) '_PS1/n' num2str(order) '_PS1.mat']);
else
    load([data '/n'  num2str(order) '_PS2/n' num2str(order) '_PS2.mat']);
end

UMean(:,n)        = UMean0_X(:,2)/1;
UMeanSigma(:,n)   = UMeanSigma_X(:,2)/1;

RMean(:,n,:)      = [RMean0_XX(:,2),RMean0_YY(:,2),RMean0_ZZ(:,2),RMean0_XY(:,2)]/ut2;
RMeanSigma(:,n,:) = [RMeanSigma_XX(:,2),RMeanSigma_YY(:,2),RMeanSigma_ZZ(:,2),RMeanSigma_XY(:,2)]/ut2;

end

caseName = 'n4_FSvsPS_comp';
uqPlotsDir =[data '/' caseName];
% uqPlotsDir = './uqPlots';

if ~exist(uqPlotsDir, 'dir')
    mkdir(uqPlotsDir);
end
if ~exist([uqPlotsDir '/tex'])
    mkdir([uqPlotsDir '/tex']);
end
if ~exist([uqPlotsDir '/pdf'])
    mkdir([uqPlotsDir '/pdf']);
end
if ~exist([uqPlotsDir '/png'])
    mkdir([uqPlotsDir '/png']);
end

% load([data '/' caseName '/' caseName '.mat']);

%% Plot settings

plotU    = 0;
plotUrms = 0;
plotuv   = 0;

plotU    = 1;
plotUrms = 1;
plotuv   = 1;

% Defaults 
uqClr  = {rgb('DarkBlue'),rgb('Navy'),rgb('Blue'),...
            rgb('RoyalBlue'),rgb('MediumSlateBlue')};
uqLS   = {'-','-','--','-.',':'};
uqMark = {'none','.','none','none','none'};
DNSclr   = '-k';
MCclr   = '-r';
LW1      = 0.5;
LW0_5    = 0.5;

fill_color = [.5 .5 .5];
FaceAlpha = 0.5;
EdgeColor = 'none';

fSize = 25;
txtFSize = 20;
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

allLeg  = {'MC','FS','PS1','PS2'};
legBool = 'off';

%% U

for makePlot = 0 : plotU-1
    
figure(1);

    ut      = 1.7*nu/yWall 
    % ut		= YPlus*nu/yWall 		% LES
    ut2		= ut*ut;
    ut_nu 	= ut/nu;
    yPlus   = y*ut_nu;

    MC_ut     = 0.0071;
    MC_factor = MC_ut/ut;

    mean = UMean/ut;
    MCmean  = UMean0_MC'*MC_factor;

    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$\mu_u /u_\tau$");
    
    plt1 = plot(MC_y, MCmean, MCclr,'LineWidth',LW0_5);

    for n = 1:numCases
        plt(n) = plot(yByd, mean(:,n), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0 1]);
    ylim([0 25]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

%     set(figure(1), 'Position',[105   647   560   420]);
%     MagInset(figure(1), -1, [0.2 0.4 15 20], [0.5 0.95 2.5 15], {'NE','NE';'SW','SW'});
    
    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','southeast');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;
%     uqPlotsPrint(uqPlotsDir,'1UMean_vs_y');
  

figure(2);

%     ut      = 1.7*nu/yWall 
    ut		= YPlus*nu/yWall 		% LES
    ut2		= ut*ut;
    ut_nu 	= ut/nu;
    yPlus   = y*ut_nu;
    
    MC_ut     = 0.0079;
    MC_factor = MC_ut/ut;
    
    sigma = UMeanSigma/ut;
    var = sigma.^2;
    MCsigma = MC_factor*UMeanSigma_MC';
    MCvar = MCsigma.^2;

    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$\sigma_u^2 /u_\tau ^2$");
    
    plt1 = plot(MC_y, MCvar, MCclr,'LineWidth',LW0_5);

    for n = 1:numCases
        plt(n) = plot(yByd, var(:,n), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0 1]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

%     ylim([0 0.4]);
%     set(figure(2), 'Position',[105   647   560   420]);
%     MagInset(figure(2), -1, [0.2 0.4 15 20], [0.5 0.95 2.5 15], {'NE','NE';'SW','SW'});

    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'1UVar_vs_y');
    
    
figure(3);

    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$(\sigma/\mu)_u$");
    
    start = 1;
    shift = 0.04;
    plt1 = plot(MC_y(start:end)-shift, MCsigma(start:end)./MCmean(start:end), MCclr,'LineWidth',LW0_5);

    for n = 1:numCases
        plt(n) = plot(yByd(start:end)-shift, sigma(start:end,n)./mean(start:end,n), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0.0 1]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

%     ylim([0 0.05]);
    
    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'1USigmaByMean_vs_y');
     
end

%% uu

for makePlot = 0 : plotUrms-1
    
mean   = sqrt(RMean);
MCmean = sqrt([RMean0_XX_MC', RMean0_YY_MC', RMean0_ZZ_MC']);

sigma     = 0.75*(RMeanSigma);
MCsigma   = ([RMeanSigma_XX_MC', RMeanSigma_YY_MC', RMeanSigma_ZZ_MC']);

offset = 1;
ymax = 6*offset;

figure(4)
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
%     ylabel("$( \mu_{v_{rms}} , \mu_{w_{rms}} + u_\tau , \mu_{u_{rms}} + 2u_\tau ) / u_\tau$")
    ylabel("$( \mu_{v_{rms}} , \mu_{w_{rms}}, \mu_{u_{rms}} ) / u_\tau$")

    %vv
    plt1 = plot(MC_y, MCmean(:,2), MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, mean(:,n,2), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    %ww
    plt1 = plot(MC_y, MCmean(:,3)+1*offset, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, mean(:,n,3)+1*offset, 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    %uu
    plt1 = plot(MC_y, MCmean(:,1)+2*offset, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, mean(:,n,1)+2*offset, 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0 1]);
    ylim([0 ymax]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    
    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2RMean_vs_y');
    

figure(5)   
    
    offset = 0.4;
    ymax = 6*offset;
    
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
%     ylabel("$( \sigma_{v_{rms}}^2, \sigma_{w_{rms}}^2 + 0.4u_\tau^2 , \sigma_{u_{rms}}^2 + 0.8u_\tau^2 ) / u_\tau^2$")
    ylabel("$( \sigma_{v_{rms}}^2, \sigma_{w_{rms}}^2, \sigma_{u_{rms}}^2 ) / u_\tau^2$")

    %vv
    plt1 = plot(MC_y, MCsigma(:,2), MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, sigma(:,n,2), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    %ww
    plt1 = plot(MC_y, MCsigma(:,3)+1*offset, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, sigma(:,n,3)+1*offset, 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    %uu
    plt1 = plot(MC_y, MCsigma(:,1)+2*offset, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, sigma(:,n,1)+2*offset, 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0 1]);
    ylim([0 ymax]);
    yticks([0 0.4 0.8 1.20 1.60 2.00 2.40]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2RMeanSigma_vs_y');
    

figure(6)   

    sigma     = sqrt(sigma);
    MCsigma   = sqrt(MCsigma);
    
    offset = 0.8;
    ymax = 3*offset;
    
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
    ylabel("$ (\sigma/\mu)_{v_{rms}},(\sigma/\mu)_{w_{rms}},(\sigma/\mu)_{u_{rms}}$")

    %vv
    plt1 = plot(MC_y(2:end), MCsigma(2:end,2)./MCmean(2:end,2), MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd(2:end), sigma(2:end,n,2)./mean(2:end,n,2), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    %ww
    plt1 = plot(MC_y(2:end), MCsigma(2:end,3)./MCmean(2:end,3)+1*offset, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd(2:end), sigma(2:end,n,3)./mean(2:end,n,3)+1*offset, 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    %uu
    plt1 = plot(MC_y(2:end), MCsigma(2:end,1)./MCmean(2:end,1)+2*offset, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd(2:end), sigma(2:end,n,1)./mean(2:end,n,1)+2*offset, 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0 1]);
    ylim([0 ymax]);
    yticks([0 0.4 0.8 1.20 1.60 2.00 2.40]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    
    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2RSigmaByMean_vs_y');
    
end

%% uv

for makePlot = 0 : plotuv-1
    
mean   = -RMean;
MCmean = -RMean0_XY_MC';

sigma     = 0.75*(RMeanSigma);
MCsigma   = RMeanSigma_XY_MC';

figure(7)
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
    ylabel("$\mu_{-\overline{u'v'}} / u_\tau^2$")

    %uv
    plt1 = plot(MC_y, MCmean, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, mean(:,n,4), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end

    xlim([0 1]);
    ylim([0 1]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    
    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2uvRMean_vs_y');
    

figure(8)   
    
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
    ylabel("$\sigma_{-\overline{u'v'}} / u_\tau^2$")

    %uv
    plt1 = plot(MC_y, MCsigma, MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd, sigma(:,n,4), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
     
    xlim([0 1]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2uvRMeanSigma_vs_y');
    

figure(9)   
   
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
    ylabel("$(\sigma/\mu)_{-\overline{u'v'}}$")

    %uv
    plt1 = plot(MC_y(2:end), MCsigma(2:end)./MCmean(2:end), MCclr,'LineWidth',LW0_5);
    for n = 1:numCases
        plt(n) = plot(yByd(2:end), sigma(2:end,n,4)./mean(2:end,n,4), 'color', uqClr{n},'LineWidth', LW0_5,...
                      'LineStyle',uqLS{n},'marker',uqMark{n});
    end
    
    xlim([0 1]);
    ylim([0 1]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt(1) plt(2) plt(3)], allLeg ,'Location','northeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2uvRSigmaByMean_vs_y');
    
end

%% Save mat file
save([data '/' caseName '/' caseName '.mat'])
    
    
    
    
    
    
    
    
    
    
    