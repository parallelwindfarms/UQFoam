function [ ] = turbChFlow_DET_withDNS_main(data, caseName, plotp, plotU,plotUrms,plotuv)
%% Make new dirs

uqPlotsDir =[data '/DET/' caseName];

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

%% Load data

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/2turbChFlow']);
addpath([SCRIPTS '/matlabScripts']);
addpath([SCRIPTS '/matlabScripts/matlab2tikz-master/src']);
addpath([data '/DNS']);

Udns = dlmread("UMean.txt", '', 2, 0);
Rdns = dlmread("RMean.txt", '', 2, 0);

addpath("./postProcessing/collapsedFields/latesTimeDir");

files = [   
            "pMean.xy";
            "UMean_X.xy";
            "UPrime2Mean_XX.xy"; "UPrime2Mean_YY.xy"; 
            "UPrime2Mean_ZZ.xy"; "UPrime2Mean_XY.xy";
        ];
    
for idx = 1:length(files)
    load(files(idx));
end

%% Initialization

col1 = 1;
col2 = 2;

delta = 1;

% M     = dlmread("postProcessing/patchExpression_yPlus/0/bottomWall", '', 1,0);
% YPlus = 1.678871;           % M1
YPlus = 8.930050e-01;       % M2
y     = pMean(:,1);
yByd  = y/delta;
yWall = y(2);

Ub		= 0.1335;
nu		= 2e-5;


UT		= 0.0079;               % DNS @ Re_t = 395
UT2		= UT*UT;
UT_nu	= UT/nu;
Y       = Udns(:,1);
Ypl     = Udns(:,2);

%% Plot settings

% Defaults 
MeanClr  = '--r';
DNSclr   = '-k';
LW1      = 0.5;
LW0_5    = 0.5;

legBool = 'off';
leg_DNS = {'DNS','DET'};
          
fSize = 15;
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

%% pMean

for makePlot = 0 : plotp-1

% ut      = 1.8*nu/yWall        % M1
ut      = 0.9*nu/yWall          % M2
% ut		= YPlus*nu/yWall 		% LES
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

p = pMean(:,col2);

figure(1)
hold on;
grid on;
xlabel("$y/\delta$");
ylabel("$\overline p$")
% axis([0 1 -5e-4 5e-4])
plot(y, p, 'r');
hold off;
uqPlotsPrint(uqPlotsDir,'0pMean_vs_y');

figure(2)
semilogx(yPlus, p, 'r');
hold on;
grid on;
xlabel("$y^+/\delta$");
ylabel("$\overline p$")
hold off;
uqPlotsPrint(uqPlotsDir,'0pMean_vs_yPl');

end

%% UMean

for makePlot = 0 : plotU-1

% ut      = 1.8*nu/yWall          % M1
ut      = 0.91*nu/yWall         % M2
% ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

dns  = Udns(:,3);
Mean = UMean_X(:,col2)/ut;

figure(4)
    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\overline u /u_\tau$");
    
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);

    axis([0 1 0 30])
    
    leg1 = legend([plt1 plt2], leg_DNS,'Location','southeast');
    set(leg1,'Visible','on');
    hold off;

    uqPlotsPrint(uqPlotsDir,'1UMean_vs_y');
    
figure(44)
    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\overline u /u_\tau$");    
    
    plt1 = plot(Ypl(2:end), dns(2:end), DNSclr,'LineWidth',LW1);
    
    plt2 = plot(yPlus(2:end), Mean(2:end), MeanClr,'LineWidth',LW0_5);
    
    set(gca,'XScale','log')

    axis([0 400 0 30])

    leg1 = legend([plt1 plt2], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'1UMean_vs_yPl');

end

%% UU

for makePlot = 0 : plotUrms-1

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

dns  = sqrt(Rdns(:,3:5));
Mean = [(UPrime2Mean_XX(:,col2)), (UPrime2Mean_YY(:,col2)), (UPrime2Mean_ZZ(:,col2))];
Mean = sqrt(Mean)/ut;

offset = 1;
ymax = 6*offset;

figure(5)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau$")

    plt1 = plot(Y, dns(:,2), DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,2), MeanClr,'LineWidth',LW0_5);
    
    plt1 = plot(Y, dns(:,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    
    plt1 = plot(Y, dns(:,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    
    axis([0 1 0 ymax])

    leg1 = legend([plt1 plt2], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'2RMeanSigma_vs_y');
    
figure(55)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau$")

    plt1 = plot(Ypl(2:end), dns(2:end,2), DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,2), MeanClr,'LineWidth',LW0_5);
    
    plt1 = plot(Ypl(2:end), dns(2:end,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    
    plt1 = plot(Ypl(2:end), dns(2:end,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    
    set(gca,'XScale','log')
    axis([0 400 0 ymax])

    leg1 = legend([plt1 plt2], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;
    

    uqPlotsPrint(uqPlotsDir,'2RMeanSigma_vs_yPl');

end

%% uv

for makePlot = 0 : plotuv-1
    
% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

dns  = -Rdns(:,6);
Mean = -UPrime2Mean_XY(:,col2)/ut2;

figure(6)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$- \overline{u'v'} /u_\tau^2$")
      
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);

    axis([0 1 0 1])

    leg1 = legend([plt1 plt2], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'3uvRMeanSigma_vs_y');
    
figure(66)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$- \overline{u'v'} /u_\tau^2$")
      
    plt1 = plot(Ypl(2:end), dns(2:end), DNSclr,'LineWidth',LW1);
    
    plt2 = plot(yPlus(2:end), Mean(2:end), MeanClr,'LineWidth',LW0_5);
    
    set(gca,'XScale','log')
    axis([0 400 0 1])

    leg1 = legend([plt1 plt2], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'3uvRMeanSigma_vs_yPl');

end

%% Save mat file
save([uqPlotsDir '/' caseName '.mat'])


end