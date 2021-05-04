function [ ] = turbChFlow_M2_withDNS_withDET_withNIPC_main(data, caseName, ...
                    plotp,plotnut, plotU,plotUrms,plotuv)
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

%DNS
Udns = dlmread("UMean.txt", '', 2, 0);
Rdns = dlmread("RMean.txt", '', 2, 0);

%DET
load([data '/DET_M2/DETdata_M2.mat']);

%NIPC
load(['./NIPC_M2_data.mat']);

%IPC
addpath("./postProcessing/collapsedFields/latesTimeDir");

files = [   
            "nutMean0.xy"; "nutMeanSigma.xy";
            "pMeanSigma.xy"; "UMeanSigma_X.xy"; 
            "RMeanSigma_XX.xy"; "RMeanSigma_YY.xy"; "RMeanSigma_ZZ.xy"; "RMeanSigma_XY.xy";       

            "pMean0.xy"; "UMean0_X.xy"; "RMean0_XX.xy"; "RMean0_YY.xy"; "RMean0_ZZ.xy"; "RMean0_XY.xy";
            "pMean1.xy"; "UMean1_X.xy"; "RMean1_XX.xy"; "RMean1_YY.xy"; "RMean1_ZZ.xy"; "RMean1_XY.xy";
            "pMean2.xy"; "UMean2_X.xy"; "RMean2_XX.xy"; "RMean2_YY.xy"; "RMean2_ZZ.xy"; "RMean2_XY.xy";
%             "pMean3.xy"; "UMean3_X.xy"; "RMean3_XX.xy"; "RMean3_YY.xy"; "RMean3_ZZ.xy"; "RMean3_XY.xy";

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
DETclr   = '--k';

MeanClr  = '-b';
fill_color = rgb('blue');%[.5 .5 .5];
FaceAlpha = 0.2;
EdgeColor = 'none';

NIPCclr   = '-r';
NIPC_fill_color = rgb('red');%[.5 .5 .5];
NIPC_FaceAlpha = 0.2;
NIPC_EdgeColor = 'none';
NIPC_N   = 2

LW1      = 0.25;
LW0_5    = 0.25;

legBool        = 'off';
leg_DNS_NIPC_IPC = {'DNS','DET','$\mathbf{E}_{\mathrm{IPC}}$','$\mathbf{E}_{\mathrm{NIPC}}$',...
              '$\pm 2\sqrt{\mathbf{V}}_{\mathrm{IPC}}$','$\pm 2\sqrt{\mathbf{V}}_{\mathrm{NIPC}}$'};
          
fSize = 25;
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', fSize);
% set(0,'DefaultLegendFontSize',fSize);
% set(0,'DefaultTextFontSize',fSize);

% set(groot, 'units', 'inches', 'position', [0 0 8 4])
set(groot, 'defaultFigureUnits','inches');
set(groot, 'defaultFigurePosition',[2.5 1.5 8 6]);

% set(groot, 'defaultFigurePaperPositionMode', 'manual');
set(groot, 'defaultFigurePaperUnits', 'inches');
set(groot, 'defaultFigurePaperPosition', [2.5 1.5 8 6]);


%%

SigmaMean = false;

%% pMean
N_p = 2
N = N_p;

for makePlot = 0 : plotp-1

p = pMean0(:,col2);
pMeanSigma = pMeanSigma(:,col2);

figure(1)
hold on;
grid off;
xlabel("$y/\delta$");
ylabel("$\overline p$")
axis([0 1 -1e-4 1e-4])
plot(y, p, 'b');
plot(y, p + N*pMeanSigma, 'g', y, p - N*pMeanSigma, 'g');
offset = 0.00053;
% plot(NIPC_y, offset + pMean0_NIPC', 'r');
% plot(NIPC_y, offset + pMean0_NIPC' + NIPC_N*pMeanSigma_NIPC', 'y', NIPC_y, offset + pMean0_NIPC' - NIPC_N*pMeanSigma_NIPC', 'y');
hold off;
uqPlotsPrint(uqPlotsDir,'0pMean_vs_y');

end

%% nutMean
N_nut = 2%0.75
N = N_nut;

for makePlot = 0 : plotnut-1

Mean  = nutMean0(:,col2);
Sigma = nutMeanSigma(:,col2);
% Sigma = sqrt(nutMean1(:,col2).^2 + nutMean2(:,col2).^2);
% Sigma = nutMean1(:,col2);
% Sigma = nutMean2(:,col2);
Ub    = Mean+N*Sigma;
Lb    = Mean-N*Sigma;

figure(2)
    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$\nu_t$");
        
    %DET
%     plt4 = plot(yByd, det, DETclr,'LineWidth',LW0_5);
    
    % IPC
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    axis([0 1 0 2e-5])
%     leg1 = legend([plt2 fl1], leg_DNS_NIPC_IPC,'Location','southeast');
%     set(leg1,'Visible','on');
%     set(leg1,'Color','none');
%     set(leg1,'EdgeColor','none');
    hold off;
    
%     set(figure(2), 'Position',[2.5 1.5 8 6])
%     set(gcf, 'PaperPosition', [2.5 1.5 8 6]);
    uqPlotsPrint(uqPlotsDir,'0nutMean_vs_y');

end

%% UMean
N_U = 1.25
N = N_U;

for makePlot = 0 : plotU-1

ut      = 0.86*nu/yWall 
% ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;
utDET  = 0.0076;

dns  = Udns(:,3);
det  = UMean_X_DET_M2(:,col2)/utDET;

Mean  = UMean0_X(:,col2)/ut;
% Sigma = UMeanSigma_X(:,col2)/ut;
% Sigma = sqrt(UMean1_X(:,col2).^2)/ut;
Sigma = sqrt(UMean1_X(:,col2).^2+UMean2_X(:,col2).^2)/ut;
% Sigma = sqrt(UMean1_X(:,col2).^2+UMean2_X(:,col2).^2+UMean3_X(:,col2).^2)/ut;
Ub    = Mean+N*Sigma;
Lb    = Mean-N*Sigma;

if (SigmaMean == true)
    Ub   = Mean+N*USigmaMean_X(:,col2)/ut;
    Lb   = Mean-N*USigmaMean_X(:,col2)/ut;
end

NIPC_ut   = 0.0073;
NIPC_Mean = UMean0_NIPC/NIPC_ut;
NIPC_Ub   = NIPC_Mean+NIPC_N*UMeanSigma_NIPC/NIPC_ut;
NIPC_Lb   = NIPC_Mean-NIPC_N*UMeanSigma_NIPC/NIPC_ut;

figure(4)
    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$\overline u /u_\tau$");
    
    %DNS
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    %DET
    plt4 = plot(yByd, det, DETclr,'LineWidth',LW0_5);

    %IPC
    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    axis([0 1 0 25])
    
    %NIPC
    plt3 = plot(NIPC_y, NIPC_Mean, NIPCclr,'LineWidth',LW0_5);
    fl2 = fill([NIPC_y', fliplr(NIPC_y')], [NIPC_Ub', fliplr(NIPC_Lb')], ...
              NIPC_fill_color, 'FaceAlpha', NIPC_FaceAlpha,'EdgeColor', NIPC_EdgeColor);
    uistack(fl2,'bottom');

    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    leg1 = legend([plt1 plt4 plt2 plt3 fl1 fl2], leg_DNS_NIPC_IPC,'Location','southeast');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'1UMean_vs_y');
 
end

%% UU
N_R = 2
N = N_R;

for makePlot = 0 : plotUrms-1

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;
ut2DET  = (0.88*nu/yWall)^2 %0.007^2;

dns  = sqrt(Rdns(:,3:5));
det  = sqrt([(UPrime2Mean_XX_DET_M2(:,col2)), (UPrime2Mean_YY_DET_M2(:,col2)), ...
             (UPrime2Mean_ZZ_DET_M2(:,col2))]/ut2DET);
         
Mean = [(RMean0_XX(:,col2)), (RMean0_YY(:,col2)), (RMean0_ZZ(:,col2))];
% Sigma2 = [(RMeanSigma_XX(:,col2)), (RMeanSigma_YY(:,col2)), (RMeanSigma_ZZ(:,col2))];
% Sigma2 = sqrt([RMean1_XX(:,col2).^2, RMean1_YY(:,col2).^2, RMean1_ZZ(:,col2).^2]);
Sigma2 = sqrt([RMean1_XX(:,col2).^2+RMean2_XX(:,col2).^2, ...
               RMean1_YY(:,col2).^2+RMean2_YY(:,col2).^2, ...
               RMean1_ZZ(:,col2).^2+RMean2_ZZ(:,col2).^2]);
% Sigma2 = sqrt([RMean1_XX(:,col2).^2+RMean2_XX(:,col2).^2+RMean3_XX(:,col2).^2, ...
%                RMean1_YY(:,col2).^2+RMean2_YY(:,col2).^2+RMean3_YY(:,col2).^2, ...
%                RMean1_ZZ(:,col2).^2+RMean2_ZZ(:,col2).^2+RMean3_ZZ(:,col2).^2]);
Ub   = sqrt(Mean+N*Sigma2)/ut;
Lb   = sqrt(Mean-N*Sigma2)/ut;

if (SigmaMean == true)
    Ub   = sqrt(Mean+N*[(RSigmaMean_XX(:,col2)), (RSigmaMean_YY(:,col2)), (RSigmaMean_ZZ(:,col2))])/ut;
    Lb   = sqrt(Mean-N*[(RSigmaMean_XX(:,col2)), (RSigmaMean_YY(:,col2)), (RSigmaMean_ZZ(:,col2))])/ut;
end
Mean = sqrt(Mean)/ut;

NIPC_ut = 0.007012753335948415;
NIPC_Mean = RMean0_NIPC(:,1:3);
NIPC_Ub   = sqrt(NIPC_Mean+NIPC_N*RMeanSigma_NIPC(:,1:3))/NIPC_ut;
NIPC_Lb   = sqrt(NIPC_Mean-NIPC_N*RMeanSigma_NIPC(:,1:3))/NIPC_ut;
NIPC_Mean = sqrt(NIPC_Mean)/NIPC_ut;

offset = 1;
ymax = 6*offset;

figure(5)
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
%     ylabel("$( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau$")
    ylabel("$( v_{rms} , w_{rms} , u_{rms} ) / u_\tau$") % For poster
    
    %DNS
    plt1 = plot(Y, dns(:,2), DNSclr,'LineWidth',LW1);
    %DET
    plt4 = plot(yByd, det(:,2), DETclr,'LineWidth',LW0_5);
    %IPC
    plt2 = plot(yByd, Mean(:,2), MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub(:,2)', fliplr(Lb(:,2)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    %NIPC
    plt3 = plot(NIPC_y, NIPC_Mean(:,2), NIPCclr,'LineWidth',LW0_5);
    fl2 = fill([NIPC_y', fliplr(NIPC_y')], [NIPC_Ub(:,2)', fliplr(NIPC_Lb(:,2)')], ...
              NIPC_fill_color, 'FaceAlpha', NIPC_FaceAlpha,'EdgeColor', NIPC_EdgeColor);
    uistack(fl2,'bottom');
    
    plt1 = plot(Y, dns(:,3) + 1*offset, DNSclr,'LineWidth',LW1);
    plt4 = plot(yByd, det(:,3) + 1*offset, DETclr,'LineWidth',LW0_5);
    plt2 = plot(yByd, Mean(:,3) + 1*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [(Ub(:,3) + 1*offset)', fliplr((Lb(:,3) + 1*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    plt3 = plot(NIPC_y, NIPC_Mean(:,3) + 1*offset, NIPCclr,'LineWidth',LW0_5);
    fl2 = fill([NIPC_y', fliplr(NIPC_y')], [(NIPC_Ub(:,3) + 1*offset)', fliplr((NIPC_Lb(:,3) + 1*offset)')], ...
              NIPC_fill_color, 'FaceAlpha', NIPC_FaceAlpha,'EdgeColor', NIPC_EdgeColor);
    uistack(fl2,'bottom');
    
    plt1 = plot(Y, dns(:,1) + 2*offset, DNSclr,'LineWidth',LW1);
    plt4 = plot(yByd, det(:,1) + 2*offset, DETclr,'LineWidth',LW0_5);
    plt2 = plot(yByd, Mean(:,1) + 2*offset, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [(Ub(:,1) + 2*offset)', fliplr((Lb(:,1) + 2*offset)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');   
    plt3 = plot(NIPC_y, NIPC_Mean(:,1) + 2*offset, NIPCclr,'LineWidth',LW0_5);
    fl2 = fill([NIPC_y', fliplr(NIPC_y')], [(NIPC_Ub(:,1) + 2*offset)', fliplr((NIPC_Lb(:,1) + 2*offset)')], ...
              NIPC_fill_color, 'FaceAlpha', NIPC_FaceAlpha,'EdgeColor', NIPC_EdgeColor);
    uistack(fl2,'bottom');
    
    axis([0 1 0 ymax-1])
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    leg1 = legend([plt1 plt4 plt2 plt3 fl1 fl2], leg_DNS_NIPC_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'2RMeanSigma_vs_y');
    
end

%% uv

for makePlot = 0 : plotuv-1
    
ut      = 0.85*nu/yWall 
% ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;
ut2DET  = (0.9*nu/yWall)^2 %0.007^2;

dns  = -Rdns(:,6);
det  = -UPrime2Mean_XY_DET_M2(:,col2)/ut2DET;

Mean  = -RMean0_XY(:,col2)/ut2;
% Sigma = RMeanSigma_XY(:,col2)/ut2;
% Sigma = sqrt(RMean1_XY(:,col2).^2)/ut2;
Sigma = sqrt(RMean1_XY(:,col2).^2+RMean2_XY(:,col2).^2)/ut2;
% Sigma = sqrt(RMean1_XY(:,col2).^2+RMean2_XY(:,col2).^2+RMean3_XY(:,col2).^2)/ut2;
Ub = Mean+N*Sigma;
Lb = Mean-N*Sigma;

if (SigmaMean == true)
    Ub   = Mean + N*RSigmaMean_XY(:,col2)/ut2;
    Lb   = Mean - N*RSigmaMean_XY(:,col2)/ut2;
end 

NIPC_ut2  = NIPC_ut^2;
NIPC_Mean = -RMean0_NIPC(:,4)/NIPC_ut2;
NIPC_Ub   = NIPC_Mean+NIPC_N*RMeanSigma_NIPC(:,4)/NIPC_ut2;
NIPC_Lb   = NIPC_Mean-NIPC_N*RMeanSigma_NIPC(:,4)/NIPC_ut2;

figure(6)
    hold on;
    grid off;
    
    xlabel("$y/\delta$");
    ylabel("$- \overline{u'v'} /u_\tau^2$")
      
    plt1 = plot(Y, dns, DNSclr,'LineWidth',LW1);
    
    plt4 = plot(yByd, det, DETclr,'LineWidth',LW0_5);

    plt2 = plot(yByd, Mean, MeanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [Ub', fliplr(Lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');

    plt3 = plot(NIPC_y, NIPC_Mean, NIPCclr,'LineWidth',LW0_5);
    fl2 = fill([NIPC_y', fliplr(NIPC_y')], [NIPC_Ub', fliplr(NIPC_Lb')], ...
              NIPC_fill_color, 'FaceAlpha', NIPC_FaceAlpha,'EdgeColor', NIPC_EdgeColor);
    uistack(fl2,'bottom');
    
    axis([0 1 0 1])    
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    leg1 = legend([plt1 plt4 plt2 plt3 fl1 fl2], leg_DNS_NIPC_IPC,'Location','northeast');
    set(leg1,'Visible','off');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'3uvRMeanSigma_vs_y');
    
end

%% Save mat file
save([uqPlotsDir '/' caseName '.mat'])


end