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

order       = 2;
% numCases    = 3;
numCases    = 2;

allCases    = {[data '/n'  num2str(order) '_FS/n'     num2str(order) '_FS.mat'],...
%                [data '/n'  num2str(order) '_M2_FS/n'  num2str(order) '_M2_FS.mat'],...
               [data '/n'  num2str(order) '_M2_PS1/n' num2str(order) '_M2_PS1.mat']};%,...
%                [data '/n'  num2str(order) '_M2_PS2/n' num2str(order) '_M2_PS2.mat']};

for i = 1:numCases

load(allCases{i});

UMean{i}      = UMean0_X(:,2)/1;
UMeanSigma{i} = UMeanSigma_X(:,2)/1;

RMean{i}      = [RMean0_XX(:,2),RMean0_YY(:,2),RMean0_ZZ(:,2),RMean0_XY(:,2)]/ut2;
RMeanSigma{i} = [RMeanSigma_XX(:,2),RMeanSigma_YY(:,2),RMeanSigma_ZZ(:,2),RMeanSigma_XY(:,2)]/ut2;

YByd{i}       = yByd;
yW{i}         = yWall;
yPl{i}        = YPlus;
Ut{i}         = ut;

end

caseName = 'n2_FSM1vsFSPSM2_comp';
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
DNSclr   = '-k';

% MeanClr  = {['- ' rgb('blue')],['- ' rgb('gray')]};
MeanClr  = {'-b','-m','--k'};
fill_color = {'b','m',rgb('gray')};
FaceAlpha = {0.2,0.2,0.2};
EdgeColor = {'none','none','none'};

LW1      = 0.5;
LW0_5    = 0.5;

legBool        = 'off';

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

allLeg  = {'DNS','$\mu_{\textrm{(M1)}}$','$\mu_{\textrm{(M2)}}$',...
                 '$\pm 2\sigma_{\textrm{(M1)}}$','$\pm 2\sigma_{\textrm{(M2)}}$'};
legBool = 'off';

%% U

for makePlot = 0 : plotU-1
    
dns  = Udns(:,3);

figure(1);

    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$\mu_u /u_\tau$");   
    
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);

    for i = 1:numCases
        if i==1
            ut  = 0.0075  
        else
            ut  = 0.0077
        end
%         ut		= yPl{i}*nu/yW{i} 		% LES
        yByd = YByd{i};
        Mean = UMean{i}/ut;
        Ub   = Mean+N*UMeanSigma{i}/ut;
        Lb   = Mean-N*UMeanSigma{i}/ut;
        plt(i)  = plot(yByd, Mean, MeanClr{i},'LineWidth', LW0_5);
        fl(i) = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
                  fill_color{i}, 'FaceAlpha', FaceAlpha{i},'EdgeColor', EdgeColor{i});
        uistack(fl(i),'bottom');
    end
    
    xlim([0 1]);
    ylim([0 30]);  
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt(1) plt(2) fl(1) fl(2)], allLeg ,'Location','southeast');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;
    uqPlotsPrint(uqPlotsDir,'1UMean_vs_y');
       
end

%% UU

N = 1.5

for makePlot = 0 : plotUrms-1

dns  = sqrt(Rdns(:,3:5));

offset = 1;
ymax = 6*offset;

figure(2);

    hold on;
    grid off;

    xlabel("$y/\delta$");
%     ylabel("$( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau$")
    ylabel("$( v_{rms} , w_{rms} , u_{rms} ) / u_\tau$")
    
    plt1 = plot(Y, dns(:,2), DNSclr,'LineWidth',LW1);
    plt1 = plot(Y, dns(:,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt1 = plot(Y, dns(:,1) + 2*offset, DNSclr,'LineWidth',LW1);


    for i = 1:numCases
        yByd = YByd{i};
        Mean = RMean{i}(:,1:3);
        Ub   = sqrt(Mean+N*(RMeanSigma{i}(:,1:3)));
        Lb   = sqrt(Mean-N*(RMeanSigma{i}(:,1:3)));
        Mean = sqrt(Mean);
        
        plt(i)  = plot(yByd, Mean(:,2), MeanClr{i},'LineWidth', LW0_5);
        fl(i) = fill([yByd', fliplr(yByd')], [Ub(:,2)', fliplr(Lb(:,2)')], ...
                  fill_color{i}, 'FaceAlpha', FaceAlpha{i},'EdgeColor', EdgeColor{i});
        uistack(fl(i),'bottom');        
        
        plt(i)  = plot(yByd, Mean(:,3)+ 1*offset, MeanClr{i},'LineWidth', LW0_5);
        fl(i) = fill([yByd', fliplr(yByd')], [(Ub(:,3)+ 1*offset)', fliplr((Lb(:,3)+ 1*offset)')], ...
                  fill_color{i}, 'FaceAlpha', FaceAlpha{i},'EdgeColor', EdgeColor{i});
        uistack(fl(i),'bottom');        
        
        plt(i)  = plot(yByd, Mean(:,1)+ 2*offset, MeanClr{i},'LineWidth', LW0_5);
        fl(i) = fill([yByd', fliplr(yByd')], [(Ub(:,1)+ 2*offset)', fliplr((Lb(:,1)+ 2*offset)')], ...
                  fill_color{i}, 'FaceAlpha', FaceAlpha{i},'EdgeColor', EdgeColor{i});
        uistack(fl(i),'bottom');
    end
    
    xlim([0 1]);
    ylim([0 6]);    
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt(1) plt(2) fl(1) fl(2)], allLeg ,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'2RMean_vs_y');
    
end

%% UV

N = 1.5

for makePlot = 0 : plotUrms-1

dns  = -Rdns(:,6);

figure(3);

    hold on;
    grid off;

    xlabel("$y/\delta$");
    ylabel("$- \overline{u'v'} /u_\tau^2$")
    
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);

    for i = 1:numCases
        yByd = YByd{i};
        Mean = -RMean{i}(:,4);
        Ub   = Mean+N*(RMeanSigma{i}(:,4));
        Lb   = Mean-N*(RMeanSigma{i}(:,4));
        
        plt(i)  = plot(yByd, Mean, MeanClr{i},'LineWidth', LW0_5);
        fl(i) = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
                  fill_color{i}, 'FaceAlpha', FaceAlpha{i},'EdgeColor', EdgeColor{i});
        uistack(fl(i),'bottom');        
    end
    
    xlim([0 1]);
    ylim([0 1]);    
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt(1) plt(2) fl(1) fl(2)], allLeg ,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;
    uqPlotsPrint(uqPlotsDir,'3uvRMean_vs_y');
    
end

%% Save mat file
save([data '/' caseName '/' caseName '.mat'])
    
    
    
    
    
    
    
    
    
    
    