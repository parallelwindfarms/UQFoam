function [ ] = turbChFlow_withDNS_main(data, caseName, plotp, plotU,plotUrms,plotuv)
%% Make new dirs

uqPlotsDir =[data '/' caseName];

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
            "pMean0.xy"; "pMeanSigma.xy";
            "UMean0_X.xy"; "UMeanSigma_X.xy";
            "RMean0_XX.xy"; "RMean0_YY.xy"; "RMean0_ZZ.xy"; "RMean0_XY.xy";
            "RMeanSigma_XX.xy"; "RMeanSigma_YY.xy"; "RMeanSigma_ZZ.xy"; "RMeanSigma_XY.xy";       
        ];
    
for idx = 1:length(files)
    load(files(idx));
end

%% Initialization

col1 = 1;
col2 = 2;

delta = 1;

M     = dlmread("postProcessing/patchExpression_yPlus/0/bottomWall", '', 1,0);
YPlus = mean(M(:,col2))
y     = pMean0(:,1);
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
MeanClr  = '-b';
fill_color = rgb('blue');%[.5 .5 .5];
FaceAlpha = 0.2;
EdgeColor = 'none';

DNSclr   = '-k';
LW1      = 0.5;
LW0_5    = 0.5;

legBool = 'off';
leg_DNS = {'DNS','$\mu_{\mathrm{IPC}}$','$\pm 2\sigma_{\mathrm{IPC}}$'};
          
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

N		= 2

for makePlot = 0 : plotp-1

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

p = pMean0(:,col2);
% pSigmaMean = pSigmaMean(:,col2);
pMeanSigma = pMeanSigma(:,col2);

figure(1)
hold on;
grid on;
xlabel("$y/\delta$");
ylabel("$\overline p$")
% axis([0 1 -5e-4 5e-4])
plot(y, p, 'r');
% plot(y, p + N*pSigmaMean, 'b', y, p - N*pSigmaMean, 'b');
plot(y, p + N*pMeanSigma, 'g', y, p - N*pMeanSigma, 'g');
% legend('U','\pm sigmaMean','','\pm meanSigma','')
hold off;
uqPlotsPrint(uqPlotsDir,'0pMean_vs_y');

figure(2)
semilogx(yPlus, p, 'r');
hold on;
grid on;
xlabel("$y^+/\delta$");
ylabel("$\overline p$")
% axis([0 400 -5e-4 5e-4])
% semilogx(yPlus, p + N*pSigmaMean, 'b', yPlus, p - N*pSigmaMean, 'b');
semilogx(yPlus, p + N*pMeanSigma, 'g', yPlus, p - N*pMeanSigma, 'g');
% legend('U','\pm sigmaMean','','\pm meanSigma','')
hold off;
uqPlotsPrint(uqPlotsDir,'0pMean_vs_yPl');

end

%% UMean

for makePlot = 0 : plotU-1

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

dns  = Udns(:,3);
Mean = UMean0_X(:,col2)/ut;
Ub   = Mean+N*UMeanSigma_X(:,col2)/ut;
Lb   = Mean-N*UMeanSigma_X(:,col2)/ut;

figure(4)
    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\overline u /u_\tau$");
    
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    axis([0 1 0 30])

    leg1 = legend([plt1 plt2 fl1], leg_DNS,'Location','southeast');
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
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [Ub(2:end)', fliplr(Lb(2:end)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    set(gca,'XScale','log')

    axis([0 400 0 30])

    leg1 = legend([plt1 plt2 fl1], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'1UMean_vs_yPl');

end

%% UU

N = 1.75

for makePlot = 0 : plotUrms-1

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

dns  = sqrt(Rdns(:,3:5));
Mean = [(RMean0_XX(:,col2)), (RMean0_YY(:,col2)), (RMean0_ZZ(:,col2))];
Ub   = sqrt(Mean+N*[(RMeanSigma_XX(:,col2)), (RMeanSigma_YY(:,col2)), (RMeanSigma_ZZ(:,col2))])/ut;
Lb   = sqrt(Mean-N*[(RMeanSigma_XX(:,col2)), (RMeanSigma_YY(:,col2)), (RMeanSigma_ZZ(:,col2))])/ut;
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
    fl1 = fill([yByd', fliplr(yByd')], [Ub(:,2)', fliplr(Lb(:,2)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    plt1 = plot(Y, dns(:,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [(Ub(:,3) + 1*offset)', fliplr((Lb(:,3) + 1*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    plt1 = plot(Y, dns(:,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [(Ub(:,1) + 2*offset)', fliplr((Lb(:,1) + 2*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');   
    
    axis([0 1 0 ymax])

    leg1 = legend([plt1 plt2 fl1], leg_DNS,'Location','southeast');
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
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [Ub(2:end,2)', fliplr(Lb(2:end,2)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    plt1 = plot(Ypl(2:end), dns(2:end,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [(Ub(2:end,3) + 1*offset)', fliplr((Lb(2:end,3) + 1*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    plt1 = plot(Ypl(2:end), dns(2:end,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [(Ub(2:end,1) + 2*offset)', fliplr((Lb(2:end,1) + 2*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');   
    
    set(gca,'XScale','log')
    axis([0 400 0 ymax])

    leg1 = legend([plt1 plt2 fl1], leg_DNS,'Location','southeast');
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
Mean = -RMean0_XY(:,col2)/ut2;
Ub   = Mean + N*RMeanSigma_XY(:,col2)/ut2;
Lb   = Mean - N*RMeanSigma_XY(:,col2)/ut2;

figure(6)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$- \overline{u'v'} /u_\tau^2$")
      
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    axis([0 1 0 1])

    leg1 = legend([plt1 plt2 fl1], leg_DNS,'Location','southeast');
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
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [Ub(2:end)', fliplr(Lb(2:end)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    set(gca,'XScale','log')
    axis([0 400 0 1])

    leg1 = legend([plt1 plt2 fl1], leg_DNS,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'3uvRMeanSigma_vs_yPl');

end

%% Save mat file
save([uqPlotsDir '/' caseName '.mat'])


end