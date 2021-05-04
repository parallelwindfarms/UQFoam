%% Clean-up

clear all;
clc;
close all;

data = getenv('turbChFlowData');

%% Load data

addpath("/home/jigar/Cartesius/cHome/scripts/matlabScripts");
addpath("/home/jigar/Cartesius/cHome/scripts/matlabScripts/matlab2tikz-master/src");

UMean       = zeros(52,5);
UMeanSigma  = zeros(52,5);
RMean       = zeros(52,5,4);
RMeanSigma  = zeros(52,5,4);

for n = 1:5

load([data '/n'  num2str(n) '_FS/n' num2str(n) '_FS.mat']);

UMean(:,n)        = UMean0_X(:,2)/1;
UMeanSigma(:,n)   = UMeanSigma_X(:,2)/1;

RMean(:,n,1)      = RMean0_XX(:,2)/ut2;
RMeanSigma(:,n,1) = RMeanSigma_XX(:,2)/ut2;

RMean(:,n,2)      = RMean0_YY(:,2)/ut2;
RMeanSigma(:,n,2) = RMeanSigma_YY(:,2)/ut2;

RMean(:,n,3)      = RMean0_ZZ(:,2)/ut2;
RMeanSigma(:,n,3) = RMeanSigma_ZZ(:,2)/ut2;

RMean(:,n,4)      = RMean0_XY(:,2)/ut2;
RMeanSigma(:,n,4) = RMeanSigma_XY(:,2)/ut2;

end

caseName = 'n12345_FScomp';
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

load([data '/' caseName '/' caseName '.mat']);

%% Plot settings

% Defaults 
meanClr  = {rgb('DarkBlue'),rgb('Navy'),...
            rgb('Blue'),rgb('RoyalBlue'),rgb('MediumSlateBlue')};
meanLS   = {'-','--','-.',':','-'};
% meanMark = {'none','.','*','o','+'};
DNSclr   = '-k';
MCclr   = '-k';
LW1      = 0.5;
LW0_5    = 0.5;

fill_color = [.5 .5 .5];
FaceAlpha = 0.5;
EdgeColor = 'none';

fSize = 15;
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', fSize);
set(0,'DefaultLegendAutoUpdate','off')
% set(0,'DefaultLegendFontSize',fSize);
% set(0,'DefaultTextFontSize',fSize);

allLeg = {'MC','$n1$','$n2$','$n3$','$n4$','$n5$'};


%% U

figure(1);

    ut      = 1.7*nu/yWall 
    % ut		= YPlus*nu/yWall 		% LES
    ut2		= ut*ut;
    ut_nu 	= ut/nu;
    yPlus   = y*ut_nu;

    mean = UMean/ut;


    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\mu_u /u_\tau$");
    
    plt1 = plot(UMean_MC(:,1), UMean_MC(:,2), MCclr,'LineWidth',LW1);

    for n = 1:5
        plt(n) = plot(yByd, mean(:,n), 'color', meanClr{n},...
                      'LineWidth', LW0_5,'LineStyle',meanLS{n});
    end
    
    xlim([0 1]);
    ylim([0 25]);
    set(figure(1), 'Position',[105   647   560   420]);
    MagInset(figure(1), -1, [0.2 0.4 15 20], [0.5 0.95 2.5 15], {'NE','NE';'SW','SW'});
    legend([plt1 plt(1) plt(2) plt(3) plt(4) plt(5)], allLeg ,'Location','Best','visible','off');
    
    hold off;
%     uqPlotsPrint(uqPlotsDir,'1UMean_vs_y');
  

figure(2);

    ut      = 1.9*nu/yWall 
    % ut		= YPlus*nu/yWall 		% LES
    ut2		= ut*ut;
    ut_nu 	= ut/nu;
    yPlus   = y*ut_nu;

    dns  = Udns(:,3);
    sigma = UMeanSigma/ut;
    var = sigma.^2;

    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\sigma_u /u_\tau$");
    
    plt1 = plot(UMeanSigma_MC(:,1), UMeanSigma_MC(:,2)*1.5, MCclr,'LineWidth',LW1);

    for n = 1:5
        plt(n) = plot(yByd, sigma(:,n), 'color', meanClr{n},...
                      'LineWidth', LW0_5,'LineStyle',meanLS{n});
    end
    
    xlim([0 1]);
    ylim([0 0.5]);
%     set(figure(2), 'Position',[105   647   560   420]);
%     MagInset(figure(2), -1, [0.2 0.4 15 20], [0.5 0.95 2.5 15], {'NE','NE';'SW','SW'});
    legend([plt1 plt(1) plt(2) plt(3) plt(4) plt(5)], allLeg ,'visible','off');
    
    hold off;
    uqPlotsPrint(uqPlotsDir,'1UVar_vs_y');
    
    
    
    
    
    
%% Save mat file
save([data '/' caseName '/' caseName '.mat'])
    
    
    
    
    
    
    
    
    
    
    