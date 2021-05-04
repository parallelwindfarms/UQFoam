function [ ] = turbChFlow_withDNS_withMC_main(data, caseName, plotp, plotU,plotUrms,plotuv)
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
addpath([SCRIPTS '/matlabScripts']);
addpath([SCRIPTS '/matlabScripts/matlab2tikz-master/src']);
addpath([data '/DNS']);

Udns = dlmread("UMean.txt", '', 2, 0);
Rdns = dlmread("RMean.txt", '', 2, 0);

load([data '/MC_data.mat']);

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
DNSclr   = '-k';

MeanClr  = '-b';
fill_color = rgb('blue');%[.5 .5 .5];
FaceAlpha = 0.2;
EdgeColor = 'none';

MCclr   = '-r';
MC_fill_color = rgb('red');%[.5 .5 .5];
MC_FaceAlpha = 0.2;
MC_EdgeColor = 'none';
MC_N      = 2

LW1      = 0.5;
LW0_5    = 0.5;

legBool        = 'off';
leg_DNS_MC_IPC = {'DNS','$\mu_{\mathrm{IPC}}$','$\mu_{\mathrm{MC}}$',...
              '$\pm 2\sigma_{\mathrm{IPC}}$','$\pm 2\sigma_{\mathrm{MC}}$'};

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
N		= 2%1.75                    % times sigma

for makePlot = 0 : plotp-1

p = pMean0(:,col2);
pMeanSigma = pMeanSigma(:,col2);

figure(1)
hold on;
grid on;
xlabel("$y/\delta$");
ylabel("$\overline p$")
axis([0 1 -1e-4 10e-4])
plot(y, p, 'b');
plot(y, p + N*pMeanSigma, 'g', y, p - N*pMeanSigma, 'g');
offset = 0.00053;
plot(MC_y, offset + pMean0_MC', 'r');
plot(MC_y, offset + pMean0_MC' + MC_N*pMeanSigma_MC', 'y', MC_y, offset + pMean0_MC' - MC_N*pMeanSigma_MC', 'y');
hold off;
uqPlotsPrint(uqPlotsDir,'0pMean_vs_y');

end

%% UMean

for makePlot = 0 : plotU-1

ut      = 1.8*nu/yWall 
% ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

MC_ut     = 0.0071;
MC_factor = MC_ut/ut;

dns  = Udns(:,3);

Mean = UMean0_X(:,col2)/ut;
Ub   = Mean+N*UMeanSigma_X(:,col2)/ut;
Lb   = Mean-N*UMeanSigma_X(:,col2)/ut;

MC_Mean = UMean0_MC'*MC_factor;
MC_Ub   = MC_Mean+MC_N*MC_factor*UMeanSigma_MC';
MC_Lb   = MC_Mean-MC_N*MC_factor*UMeanSigma_MC';

figure(4)
    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\overline u /u_\tau$");
    
    % DNS
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    % IPC
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    %MC
    plt3 = plot(MC_y, MC_Mean, MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_y', fliplr(MC_y')], [MC_Ub', fliplr(MC_Lb')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    axis([0 1 0 30])
    leg1 = legend([plt1 plt2 plt3 fl1 fl2], leg_DNS_MC_IPC,'Location','southeast');
    set(leg1,'Visible','on');
    hold off;

    uqPlotsPrint(uqPlotsDir,'1UMean_vs_y');
    
figure(44)
    hold on;
    grid on;

    xlabel("$y/\delta$"');
    ylabel("$\overline u /u_\tau$");    
    
    %DNS
    plt1 = plot(Ypl(2:end), dns(2:end), DNSclr,'LineWidth',LW1);
    
    %IPC
    plt2 = plot(yPlus(2:end), Mean(2:end), MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [Ub(2:end)', fliplr(Lb(2:end)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    set(gca,'XScale','log')

    %MC
    plt3 = plot(MC_yPlus(2:end), MC_Mean(2:end), MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_yPlus(2:end)', fliplr(MC_yPlus(2:end)')], [MC_Ub(2:end)', fliplr(MC_Lb(2:end)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    set(gca,'XScale','log')
    
    axis([0 400 0 30])
    leg1 = legend([plt1 plt2 plt3 fl1 fl2], leg_DNS_MC_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'1UMean_vs_yPl');

end

%% UU

N = 1.5

for makePlot = 0 : plotUrms-1

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

% MC_ut     = 0.0071;
% MC_factor = MC_ut/ut;

dns  = sqrt(Rdns(:,3:5));

Mean = [(RMean0_XX(:,col2)), (RMean0_YY(:,col2)), (RMean0_ZZ(:,col2))];
Ub   = sqrt(Mean+N*[(RMeanSigma_XX(:,col2)), (RMeanSigma_YY(:,col2)), (RMeanSigma_ZZ(:,col2))])/ut;
Lb   = sqrt(Mean-N*[(RMeanSigma_XX(:,col2)), (RMeanSigma_YY(:,col2)), (RMeanSigma_ZZ(:,col2))])/ut;
Mean = sqrt(Mean)/ut;

MC_Mean = [RMean0_XX_MC', RMean0_YY_MC', RMean0_ZZ_MC'];
MC_Ub   = sqrt(MC_Mean+MC_N*[RMeanSigma_XX_MC', RMeanSigma_YY_MC', RMeanSigma_ZZ_MC']);
MC_Lb   = sqrt(MC_Mean-MC_N*[RMeanSigma_XX_MC', RMeanSigma_YY_MC', RMeanSigma_ZZ_MC']);
MC_Mean = sqrt(MC_Mean);

offset = 1;
ymax = 6*offset;

figure(5)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau$")

    %vv
    plt1 = plot(Y, dns(:,2), DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,2), MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub(:,2)', fliplr(Lb(:,2)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    plt3 = plot(MC_y, MC_Mean(:,2), MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_y', fliplr(MC_y')], [MC_Ub(:,2)', fliplr(MC_Lb(:,2)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    %ww
    plt1 = plot(Y, dns(:,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [(Ub(:,3) + 1*offset)', fliplr((Lb(:,3) + 1*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    plt3 = plot(MC_y, MC_Mean(:,3) + 1*offset, MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_y', fliplr(MC_y')], [(MC_Ub(:,3) + 1*offset)', fliplr((MC_Lb(:,3) + 1*offset)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    %uu
    plt1 = plot(Y, dns(:,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yByd, Mean(:,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [(Ub(:,1) + 2*offset)', fliplr((Lb(:,1) + 2*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');  
    plt3 = plot(MC_y, MC_Mean(:,1) + 2*offset, MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_y', fliplr(MC_y')], [(MC_Ub(:,1) + 2*offset)', fliplr((MC_Lb(:,1) + 2*offset)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    axis([0 1 0 ymax])
%     legend([plt1 plt2 fl1],{'DNS','$\mu_{u_{rms}}$','$\pm 2\sigma_{u_{rms}}$'});
    leg1 = legend([plt1 plt2 plt3 fl1 fl2], leg_DNS_MC_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'2RMeanSigma_vs_y');
    
figure(55)
    hold on;
    grid on;
    
    xlabel("$y/\delta$");
    ylabel("$( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau$")

    %vv
    plt1 = plot(Ypl(2:end), dns(2:end,2), DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,2), MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [Ub(2:end,2)', fliplr(Lb(2:end,2)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    plt3 = plot(MC_yPlus(2:end), MC_Mean(2:end,2), MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_yPlus(2:end)', fliplr(MC_yPlus(2:end)')], [MC_Ub(2:end,2)', fliplr(MC_Lb(2:end,2)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    %ww
    plt1 = plot(Ypl(2:end), dns(2:end,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [(Ub(2:end,3) + 1*offset)', fliplr((Lb(2:end,3) + 1*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    plt3 = plot(MC_yPlus(2:end), MC_Mean(2:end,3) + 1*offset, MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_yPlus(2:end)', fliplr(MC_yPlus(2:end)')], [(MC_Ub(2:end,3) + 1*offset)', fliplr((MC_Lb(2:end,3) + 1*offset)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    %uu
    plt1 = plot(Ypl(2:end), dns(2:end,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt2 = plot(yPlus(2:end), Mean(2:end,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yPlus(2:end)', fliplr(yPlus(2:end)')], [(Ub(2:end,1) + 2*offset)', fliplr((Lb(2:end,1) + 2*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');   
    plt3 = plot(MC_yPlus(2:end), MC_Mean(2:end,1) + 2*offset, MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_yPlus(2:end)', fliplr(MC_yPlus(2:end)')], [(MC_Ub(2:end,1) + 2*offset)', fliplr((MC_Lb(2:end,1) + 2*offset)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    set(gca,'XScale','log')
    axis([0 400 0 ymax])
%     legend([plt1 plt2 fl1],{'DNS','$\mu_{u_{rms}}$','$\pm 2\sigma_{u_{rms}}$'});
    leg1 = legend([plt1 plt2 plt3 fl1 fl2], leg_DNS_MC_IPC,'Location','southeast');
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

MC_Mean = -RMean0_XY_MC';
MC_Ub   = MC_Mean+MC_N*RMeanSigma_XY_MC';
MC_Lb   = MC_Mean-MC_N*RMeanSigma_XY_MC';

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
    
    plt3 = plot(MC_y, MC_Mean, MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_y', fliplr(MC_y')], [MC_Ub', fliplr(MC_Lb')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    axis([0 1 0 1])
%     legend([plt1 plt2 fl1],{'DNS','$\mu_u$','$\pm 2\sigma_u$'});
    leg1 = legend([plt1 plt2 plt3 fl1 fl2], leg_DNS_MC_IPC,'Location','southeast');
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

    plt3 = plot(MC_yPlus(2:end), MC_Mean(2:end), MCclr,'LineWidth',LW0_5);
    fl2 = fill([MC_yPlus(2:end)', fliplr(MC_yPlus(2:end)')], [MC_Ub(2:end)', fliplr(MC_Lb(2:end)')], ...
              MC_fill_color, 'FaceAlpha', MC_FaceAlpha,'EdgeColor', MC_EdgeColor);
    uistack(fl2,'bottom');
    
    set(gca,'XScale','log')
    axis([0 400 0 1])
%     legend([plt1 plt2 fl1],{'DNS','$\mu_u$','$\pm 2\sigma_u$'},'Location','Best');
    leg1 = legend([plt1 plt2 plt3 fl1 fl2], leg_DNS_MC_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    hold off;

    uqPlotsPrint(uqPlotsDir,'3uvRMeanSigma_vs_yPl');

end

%% Save mat file

save([uqPlotsDir '/' caseName '.mat'])


end