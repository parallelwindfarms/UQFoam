function [ ] = flowOvrCyl_main_withDet(data, caseName, plotU_c, plotU_xByd, plotV_xByd, plotUU_xByd, plotVV_xByd, plotCd, plotCl, plotCp)
%% Make new dirs

uqPlotsDir = [data '/' caseName];

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

%% Load old data and all the scripts

SCRIPTS = getenv('SCRIPTS');
addpath([SCRIPTS '/3flowOvrCyl']);
addpath([SCRIPTS '/matlabScripts']);
addpath([SCRIPTS '/matlabScripts/matlab2tikz-master/src']);
addpath([SCRIPTS '/matlabScripts/export_fig/export_fig-master']);

% EXP data
EXPDataDir = [data '/0EXP'];
EXP_num_max = 3;
sample_xByd = {'1.06','1.54','2.02','4.00','7.00','10.0'}; 
sample = {'1_06','1_54','2_02','4_00','7_00','10_0'};  
line_num_max  = length(sample);
for EXP_num=1:EXP_num_max
    Uc_EXP{EXP_num} = load([EXPDataDir '/Uc_EXP' num2str(EXP_num) '.csv']);
        for line_num=1:line_num_max
            try
            U_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/U_xByd_EXP' ...
                                    num2str(EXP_num) '_' sample{line_num} '.csv']);
            V_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/V_xByd_EXP' ...
                                    num2str(EXP_num) '_' sample{line_num} '.csv']);
            UU_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/UU_xByd_EXP' ...
                                    num2str(EXP_num) '_' sample{line_num} '.csv']);
            VV_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/VV_xByd_EXP' ...
                                    num2str(EXP_num) '_' sample{line_num} '.csv']);
            catch
            end
        end
end
load([EXPDataDir '/Cp_EXP.csv']);

% DET data
load([data '/1detCase_M1/DETdata.mat']);

% UQ data
sampleTimeDir = './postProcessing/sample/latestTimeDir';
files = dir([sampleTimeDir '/*xByd']);
for i=1:length(files)
    load([sampleTimeDir '/' files(i).name]);
end
files = dir([sampleTimeDir '/*c']);
for i=1:length(files)
    load([sampleTimeDir '/' files(i).name]);
end
% Uc(:,1:4) = load([data '/UMeanSpAvg_smooth.xy']);
% Uc(234:966,2) = linspace(0.4018,0.440,966-234+1); % Uin
% Uc(234:966,2) = linspace(0.4035,0.440,966-234+1); % Cs
% Uc_DET(234:966,2) = linspace(0.4033,0.440,966-234+1);

load([uqPlotsDir '/forceCoeffs.mat']);
t        = coefficient200(:,1);
Cd       = coefficient200(:,3);
Cd_mean  = coefficient200_mean(:,3);
Cd_sigma = coefficient200_mean_sigma(:,3);
Cl       = coefficient200(:,4);
Cl_mean  = coefficient200_mean(:,4);
Cl_sigma = coefficient200_mean_sigma(:,4);

D           = 0.1;
R           = D/2;
x           = pc(:,1);
y           = p_xByd(:,1);
xByd        = x/D;
yByd        = y/D;

%% Initialization

rhoInf= 1.10;
Uinf  = 0.59;
Uinf2 = Uinf^2;
uCol  = 3;
RCol  = 6;
yuCol = 7;
yRCol = 13;
xDirIdx = 2;
yDirIdx = 3;

%% Plot settings

% Defaults 
meanClr = '-b';
fill_color = rgb('blue');%[.5 .5 .5];
FaceAlpha = 0.2;
EdgeColor = 'none';

tclr   = '-r';
t_fill_color = rgb('red');%[.5 .5 .5];
t_FaceAlpha = 0.4;
t_EdgeColor = 'none';

EXPclr   = {'d','^','o'};
EXPfill  = {rgb('Black'),rgb('Gray'),rgb('DarkSlateGray')};
EXPedge  = {'none','none','none'};
EXPsize  = 5;
LW2      = 2;
LW1      = 1;
LW0_5    = 0.5;

DETclr   = '--k';

legBool    = 'off';
leg_EXP_IPC = {'EXP 1','EXP 2','EXP 3', 'DET','$\mathbf{E}[\cdot]$','$\pm 2\sqrt{\mathbf{V}[\cdot]}$'};

fSize = 24;
txtFSize = 15;
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

%% U_c

%N = 3;
N = 2;

for makePlot = 0 : plotU_c-1

Mean = Uc(:,2)/Uinf;
ub   = Mean+N*Uc(:,2+uCol)/Uinf;
lb   = Mean-N*Uc(:,2+uCol)/Uinf;
xByd = xByd;

figure(1)
    hold on;
    grid off;

    xlabel("$x/D$"');
    ylabel("$u / U_{\infty}$");
        
    for EXP_num=1:EXP_num_max
        plt1(EXP_num) = plot(Uc_EXP{EXP_num}(:,1), Uc_EXP{EXP_num}(:,2), ...
                        EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
                        'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
    end 

    halfWay = length(Mean(:,1))/2;

    plt2 = plot(xByd(1:halfWay,1), Mean(1:halfWay,1), meanClr,'LineWidth',LW0_5);
    fl1  = fill([xByd(1:halfWay,1)', fliplr(xByd(1:halfWay,1)')], [ub(1:halfWay,1)', fliplr(lb(1:halfWay,1)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');

    plt2 = plot(xByd(halfWay+1:end,1), Mean(halfWay+1:end,1), meanClr,'LineWidth',LW0_5);
    fl1  = fill([xByd(halfWay+1:end,1)', fliplr(xByd(halfWay+1:end,1)')], [ub(halfWay+1:end,1)', fliplr(lb(halfWay+1:end,1)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    
    plt3 = plot(Uc_DET(:,1)/D, Uc_DET(:,2)/Uinf, DETclr,'LineWidth',LW0_5);
    
    p = calculateEllipse(0, 0, 1/2*8/8, 1/2*Uinf/2*6/8, 0);
    patch(p(:,1), p(:,2),rgb('grey'),'EdgeColor', rgb('grey'),'FaceAlpha', 0.5); 
    
    axis([-5 5 -0.4 1.2])
    xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5]);
    yticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southwest');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

%     uqPlotsPrint(uqPlotsDir,'U_centerline');
    uqPlotsPrint(uqPlotsDir,'U_centerline_smooth');
end

clear Mean; clear ub; clear lb;

%% U_xByD_1.06,1.54,2.02,4,7,10

N = 2;

for makePlot = 0 : plotU_xByd-1

plotOffset = [4.5,3.0,2.0,1.0,0.75,0.2]-1;  
textOffset = [4.75,3.2,2.2,1.15,0.85,0.28];

for line_num=1:line_num_max
    Mean(:,line_num) = U_xByd(:, xDirIdx + (line_num-1)*yuCol)/Uinf + plotOffset(line_num);
    ub(:,line_num)   = Mean(:,line_num)+N*U_xByd(:, xDirIdx + uCol + (line_num-1)*yuCol)/Uinf;
    lb(:,line_num)   = Mean(:,line_num)-N*U_xByd(:, xDirIdx + uCol + (line_num-1)*yuCol)/Uinf;
    DET(:,line_num)  = U_xByd_DET(:, xDirIdx + (line_num-1)*(uCol+1))/Uinf + plotOffset(line_num);
end

figure(2)
    hold on;
    grid off;

    xlabel("$y/D$"');
    ylabel("$u / U_{\infty}$");
    
    for line_num=1:line_num_max
        
        for EXP_num=1:EXP_num_max
            try
            plt1(EXP_num) = plot(U_xByd_EXP{EXP_num,line_num}(:,1), U_xByd_EXP{EXP_num,line_num}(:,2), ...
                            EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
                            'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
            catch
            end
        end
        
        plt2 = plot(yByd, Mean(:,line_num), meanClr,'LineWidth',LW0_5);
        fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
                  fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
        uistack(fl1,'bottom'); 
        
        plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);

        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    axis([-3 3 -1 5]);
    xticks([-3 -2 -1 0 1 2 3]);
%     yticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'U_xByd');
end

%% V_xByD_1.06,1.54,2.02,4,7,10

for makePlot = 0 : plotV_xByd-1

plotOffset = [-0.7,-1.4,-2.1,-2.5,-2.75,-2.9];  
textOffset = plotOffset + 0.075;

for line_num=1:line_num_max
    Mean(:,line_num) = U_xByd(:, yDirIdx + (line_num-1)*yuCol)/Uinf + plotOffset(line_num);
    ub(:,line_num)   = Mean(:,line_num)+N*U_xByd(:, yDirIdx + uCol + (line_num-1)*yuCol)/Uinf;
    lb(:,line_num)   = Mean(:,line_num)-N*U_xByd(:, yDirIdx + uCol + (line_num-1)*yuCol)/Uinf;
    DET(:,line_num)  = U_xByd_DET(:, yDirIdx + (line_num-1)*(uCol+1))/Uinf + plotOffset(line_num);    
end

figure(3)
    hold on;
    grid off;

    xlabel("$y/D$"');
    ylabel("$v / U_{\infty}$");
    
    for line_num=1:line_num_max
        
        for EXP_num=1:EXP_num_max
            try
            plt1(EXP_num) = plot(V_xByd_EXP{EXP_num,line_num}(:,1), V_xByd_EXP{EXP_num,line_num}(:,2), ...
                            EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
                            'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
            catch
            end
        end
        
        plt2 = plot(yByd, Mean(:,line_num), meanClr,'LineWidth',LW0_5);
        fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
                  fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
        uistack(fl1,'bottom'); 

        plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);

        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    axis([-3 3 -3 0]);
    xticks([-3 -2 -1 0 1 2 3]);
%     yticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'V_xByd');
end

%% UU_xByD_1.06,1.54,2.02,4,7,10
 
%N = [2,3,3,3,3,3];
N = [2,2,2,2,2,2];

for makePlot = 0 : plotUU_xByd-1

plotOffset = [-0.3,-0.55,-0.8,-0.95,-1.1,-1.25];  
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    Mean(:,line_num) = R_xByd(:, xDirIdx + (line_num-1)*yRCol)/Uinf2 + plotOffset(line_num);
    ub(:,line_num)   = Mean(:,line_num)+N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
    lb(:,line_num)   = Mean(:,line_num)-N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
    DET(:,line_num)  = R_xByd_DET(:, xDirIdx + (line_num-1)*(RCol+1))/Uinf2 + plotOffset(line_num);    
end

figure(4)
    hold on;
    grid off;

    xlabel("$y/D$"');
    ylabel("$u'u' / U_{\infty}^2$");
    
    for line_num=1:line_num_max
        
        for EXP_num=1:EXP_num_max
            try
            plt1(EXP_num) = plot(UU_xByd_EXP{EXP_num,line_num}(:,1), UU_xByd_EXP{EXP_num,line_num}(:,2), ...
                            EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
                            'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
            catch
            end
        end
        
        plt2 = plot(yByd, Mean(:,line_num), meanClr,'LineWidth',LW0_5);
        fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
                  fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
        uistack(fl1,'bottom'); 
        
        plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);

        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    axis([-3 3 -1.4 0]);
    xticks([-3 -2 -1 0 1 2 3]);
    yticks([-1.4 -1.2 -1.0 -0.8 -0.6 -0.4 -0.2 0]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'UU_xByd');
end

%% VV_xByD_1.06,1.54,2.02,4,7,10

%N = [2,3,3,3,3,3];
N = [2,2,2,2,2,2];
factor = [1.5,1.5,1,1,1,1];

for makePlot = 0 : plotVV_xByd-1

EXPOffset = [0.2,0,0,0,0,0];
plotOffset = [-0.4,-0.9,-1.4,-1.6,-1.75,-2.05];  
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    Mean(:,line_num) = R_xByd(:, 2*xDirIdx+1 + (line_num-1)*yRCol)/Uinf2/factor(line_num) + plotOffset(line_num);
    ub(:,line_num)   = Mean(:,line_num)+N(line_num)*R_xByd(:, 2*xDirIdx+1 + RCol + (line_num-1)*yRCol)/Uinf2;
    lb(:,line_num)   = Mean(:,line_num)-N(line_num)*R_xByd(:, 2*xDirIdx+1 + RCol + (line_num-1)*yRCol)/Uinf2;
    DET(:,line_num)  = R_xByd_DET(:, 2*xDirIdx+1 + (line_num-1)*(RCol+1))/Uinf2/factor(line_num) + plotOffset(line_num);    
end

figure(5)
    hold on;
    grid off;

    xlabel("$y/D$"');
    ylabel("$v'v' / U_{\infty}^2$");
    
    for line_num=1:line_num_max
        
        for EXP_num=1:EXP_num_max
            try
            plt1(EXP_num) = plot(VV_xByd_EXP{EXP_num,line_num}(:,1), VV_xByd_EXP{EXP_num,line_num}(:,2) + EXPOffset(line_num), ...
                            EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
                            'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
            catch
            end
        end
        
        plt2 = plot(yByd, Mean(:,line_num), meanClr,'LineWidth',LW0_5);
        fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
                  fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
        uistack(fl1,'bottom'); 
        
        plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);

        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    axis([-3 3 -2.25 0])
    xticks([-3 -2 -1 0 1 2 3]);
    yticks([-2.25 -1.75 -1.25 -0.75 -0.25]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'VV_xByd');
end

%% Cd

N = 2;

% set(0, 'defaultAxesFontSize', 36);

for makePlot = 0 : plotCd-1

TUinf_D = t*Uinf/D;

figure(6)
    hold on;
    grid off;

    xlabel("$tU_{\infty}/D$"');
    ylabel("$C_d$");
      
    Mean = Cd;
    ub   = Mean+N*Cd_sigma;
    lb   = Mean-N*Cd_sigma;

    plt1 = plot(TUinf_D, Mean, tclr,'LineWidth',LW1);
    fl1  = fill([TUinf_D', fliplr(TUinf_D')], [ub', fliplr(lb')], ...
              t_fill_color, 'FaceAlpha', t_FaceAlpha,'EdgeColor', t_EdgeColor);
    uistack(fl1,'bottom');
    
    Mean = Cd_mean;
    ub   = Mean+N*Cd_sigma;
    lb   = Mean-N*Cd_sigma;
    
    plt2 = plot(TUinf_D, Mean, meanClr,'LineWidth',LW1);
    fl2  = fill([TUinf_D', fliplr(TUinf_D')], [ub', fliplr(lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl2,'bottom');

    ylim([1 1.6]);
%     yticks([1 1.1 1.2 1.3 1.4 1.5]);
    xlim([1200 1700]);
    xticks(linspace(1200,1700,6));

    leg1 = legend([plt1 fl1 plt2 fl2],{'$\mathbf{E}[C_d]$','$\pm 2\sqrt{\mathbf{V}[C_d]}$',...
           '$\mathbf{E}[\overline{C_d}]$','$\pm 2\sqrt{\mathbf{V}[\overline{C_d}]}$'},'Location','northeast');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');   
    hold off;

    fprintf('Cd_mean       = %f\n', mean(Cd_mean(:,1)))
    fprintf('Cd_sigma      = %f\n', N*abs(mean(Cd_sigma(:,1))))
    fprintf('Cd_simga_mean = %f %s\n', 100*N*abs(mean(Cd_sigma(:,1)))/mean(Cd_mean(:,1)), '%')

%     set(figure(6), 'Units','inches')
%     set(figure(6), 'Position',[2.5 1.5 4 3])
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperPosition', [2.5 1.5 4 3]);
    uqPlotsPrint(uqPlotsDir,'Cd');
end

%% Cl

N = 2;

for makePlot = 0 : plotCl-1

TUinf_D = t*Uinf/D;

figure(7)
    hold on;
    grid off;

    xlabel("$tU_{\infty}/D$"');
    ylabel("$C_l$");
      
    Mean = Cl;
    ub   = Mean+N*Cd_sigma;
    lb   = Mean-N*Cd_sigma;

    plt1 = plot(TUinf_D, Mean, tclr,'LineWidth',LW1);
    fl1  = fill([TUinf_D', fliplr(TUinf_D')], [ub', fliplr(lb')], ...
              t_fill_color, 'FaceAlpha', t_FaceAlpha,'EdgeColor', t_EdgeColor);
    uistack(fl1,'bottom');
    
    Mean = Cl_mean;
    ub   = Mean+N*Cd_sigma;
    lb   = Mean-N*Cd_sigma;
    
    plt2 = plot(TUinf_D, Mean, meanClr,'LineWidth',LW1);
    fl2  = fill([TUinf_D', fliplr(TUinf_D')], [ub', fliplr(lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl2,'bottom');

    ylim([-1 3]);
    xlim([1200 1700]);
    xticks(linspace(1200,1700,6));

leg1 = legend([plt1 fl1 plt2 fl2],{'$\mathbf{E}[C_l]$','$\pm 2\sqrt{\mathbf{V}[C_l]}$',...
           '$\mathbf{E}[\overline{C_l}]$','$\pm 2\sqrt{\mathbf{V}[\overline{C_l}]}$'},'Location','northeast');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');   
    hold off;

    fprintf('Cl_mean       = %f\n', mean(Cl_mean(:,1)))
    fprintf('Cl_sigma      = %f\n', N*(mean(Cl_sigma(:,1))))
    fprintf('Cl_simga_mean = %f %s\n', 100*N*(mean(Cl_sigma(:,1)))/mean(Cl_mean(:,1)), '%')

%     set(figure(7), 'Units','inches')
%     set(figure(7), 'Position',[2.5 1.5 4 3])
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperPosition', [2.5 1.5 4 3]);
    uqPlotsPrint(uqPlotsDir,'Cl');
end

%% Cp

N = 2;

for makePlot = 0 : plotCp-1
  
p_mean = load([uqPlotsDir '/p_mean_cyl.csv']);
p_sigma= load([uqPlotsDir '/p_sigma_cyl.csv']);

Cp_factor = 2/(rhoInf*Uinf2);
Mean  = Cp_factor*p_mean(:,2);
ub    = Mean+N*Cp_factor*p_sigma(:,2);
lb    = Mean-N*Cp_factor*p_sigma(:,2);
theta = 180*(1-p_mean(:,1)/max(p_mean(:,1)));

figure(8)
    hold on;
    grid off;

    xlabel("$\theta$"');
    ylabel("$C_p$");

    plt1 = plot(Cp_EXP(:,1), Cp_EXP(:,2), 'ks','MarkerFaceColor',...
        EXPfill{1}, 'MarkerSize',7.5,'MarkerEdgeColor',EXPedge{1});
    
    plt2 = plot(theta, Mean, meanClr,'LineWidth',LW1);
    fl1  = fill([theta', fliplr(theta')], [ub', fliplr(lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');

    ylim([-1.5 1.5]);
    yticks([-1.5 -1 -0.5 0 0.5 1 1.5]);
    xlim([0 180]);
    xticks([0 30 60 90 120 150 180]);

    leg1 = legend([plt1 plt2 fl1],{'EXP','$\mathbf{E}[C_p]$','$\pm 2\sqrt{\mathbf{V}[C_p]}$'}, 'Location','northeast');
    set(leg1,'Visible','on');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');   
    hold off;
    
%     set(figure(8), 'Units','inches')
%     set(figure(8), 'Position',[2.5 1.5 16 6])
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperPosition', [2.5 1.5 16 6]);
    uqPlotsPrint(uqPlotsDir,'Cp');
end

%% Save mat file

save([uqPlotsDir '/' caseName '.mat'])


end