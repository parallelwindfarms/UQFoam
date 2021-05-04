function [ ] = periodicHill_main_det(data, caseName, plotU_xByd, ...
                    plotUU_xByd, plotVV_xByd, plotUV_xByd)
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
addpath([SCRIPTS '/7periodicHill']);
addpath([SCRIPTS '/matlabScripts']);
addpath([SCRIPTS '/matlabScripts/matlab2tikz-master/src']);
addpath([SCRIPTS '/matlabScripts/export_fig/export_fig-master']);

% EXP data
% EXPDataDir = [data '/0EXP'];
% EXP_num_max = 3;
sample_xByd = {'0.05','0.50','1.00','2.00','3.00',...
               '4.00','5.00','6.00','7.00','8.00'}; 
sample = {'0_05','0_50','1_00','2_00','3_00',...
          '4_00','5_00','6_00','7_00','8_00'};   
line_num_max  = length(sample);
% for EXP_num=1:EXP_num_max
%     Uc_EXP{EXP_num} = load([EXPDataDir '/Uc_EXP' num2str(EXP_num) '.csv']);
%         for line_num=1:line_num_max
%             try
%             U_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/U_xByd_EXP' ...
%                                     num2str(EXP_num) '_' sample{line_num} '.csv']);
%             V_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/V_xByd_EXP' ...
%                                     num2str(EXP_num) '_' sample{line_num} '.csv']);
%             UU_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/UU_xByd_EXP' ...
%                                     num2str(EXP_num) '_' sample{line_num} '.csv']);
%             VV_xByd_EXP{EXP_num,line_num} = load([EXPDataDir '/VV_xByd_EXP' ...
%                                     num2str(EXP_num) '_' sample{line_num} '.csv']);
%             catch
%             end
%         end
% end
% load([EXPDataDir '/Cp_EXP.csv']);

% DET data
% load([data '/1detCase_M1/DETdata.mat']);

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

H    = 28/1000;
D    = H;
R    = D/2;
x    = pc(:,1);
y    = p_xByd(:,1);
xByd = x/D;
% yByd = y/D;

% Periodic Hill Geometry
Nx = 500;
[hillX, hillY] = periodicHill_geom(H*1000, Nx);
hillX = hillX/H/1000;
hillY = hillY/H/1000;

hillClr = '-k';
hillLW = 1;

%% Initialization

rhoInf= 1.10;
Uinf  = 1;
Uinf2 = Uinf^2;
uCol  = 3;
RCol  = 6;
yuCol = 4;%7;
yRCol = 7;%13;
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
set(groot, 'defaultFigurePosition',[2.5 1.5 12 6])

% set(groot, 'defaultFigurePaperPositionMode', 'manual');
set(groot, 'defaultFigurePaperUnits', 'inches');
set(groot, 'defaultFigurePaperPosition', [2.5 1.5 12 6]);

%% U_xByD_...

N = 2;

for makePlot = 0 : plotU_xByd-1

plotOffset = [0.05, 0.5, 1,2,3,4,5,6,7,8];  
textOffset = [4.75,3.2,2.2,1.15,0.85,0.28];

for line_num=1:line_num_max
    yByd(:,line_num) = U_xByd(:, 1 + (line_num-1)*yuCol)/D;
    Mean(:,line_num) = N*U_xByd(:, xDirIdx + (line_num-1)*yuCol)/Uinf + plotOffset(line_num);
%     ub(:,line_num)   = Mean(:,line_num)+N*U_xByd(:, xDirIdx + uCol + (line_num-1)*yuCol)/Uinf;
%     lb(:,line_num)   = Mean(:,line_num)-N*U_xByd(:, xDirIdx + uCol + (line_num-1)*yuCol)/Uinf;
%     DET(:,line_num)  = U_xByd_DET(:, xDirIdx + (line_num-1)*(uCol+1))/Uinf + plotOffset(line_num);
end

figure(1)
    hold on;
    grid on;

    ylabel("$y/H$"');
    xlabel("$2u/u_{b} + x/H$");
    
    for line_num=1:line_num_max
        
%         for EXP_num=1:EXP_num_max
%             try
%             plt1(EXP_num) = plot(U_xByd_EXP{EXP_num,line_num}(:,1), U_xByd_EXP{EXP_num,line_num}(:,2), ...
%                             EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
%                             'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
%             catch
%             end
%         end
        
        plt2 = plot(Mean(:,line_num), yByd(:,line_num), meanClr,'LineWidth',LW0_5);
%         fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
%                   fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%         uistack(fl1,'bottom'); 
        
%         plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);

%         text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    pltHill = plot(hillX, hillY, hillClr,'LineWidth', hillLW);
    pltX0   = plot(0*hillX, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltx1   = plot(0*hillX+9, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltY3H  = plot(hillX, 0*hillY+85/H/1000, hillClr,'LineWidth', hillLW);
    
%     axis([-3 3 -1 5]);
    xticks([0 0.5 1 2 3 4 5 6 7 8 9 10]);
% %     yticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2]);
%     leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
%     set(leg1,'Visible',legBool);
%     set(leg1,'Color','none');
%     set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'U_xByd');
end

%% UU_xByD_...
 
%N = [2,3,3,3,3,3];
N = [2,2,2,2,2,2];

for makePlot = 0 : plotUU_xByd-1

plotOffset = [0.05, 0.5, 1,2,3,4,5,6,7,8]; 
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    yByd(:,line_num) = R_xByd(:, 1 + (line_num-1)*yRCol)/D;
    Mean(:,line_num) = 10*R_xByd(:, xDirIdx + (line_num-1)*yRCol)/Uinf2 + plotOffset(line_num);
%     ub(:,line_num)   = Mean(:,line_num)+N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
%     lb(:,line_num)   = Mean(:,line_num)-N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
%     DET(:,line_num)  = R_xByd_DET(:, xDirIdx + (line_num-1)*(RCol+1))/Uinf2 + plotOffset(line_num);    
end

figure(2)
    hold on;
    grid on;

    ylabel("$y/H$"');
    xlabel("$10u'u' + x/H$");
    
    for line_num=1:line_num_max
        
%         for EXP_num=1:EXP_num_max
%             try
%             plt1(EXP_num) = plot(UU_xByd_EXP{EXP_num,line_num}(:,1), UU_xByd_EXP{EXP_num,line_num}(:,2), ...
%                             EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
%                             'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
%             catch
%             end
%         end
        
        plt2 = plot(Mean(:,line_num),yByd(:,line_num), meanClr,'LineWidth',LW0_5);
%         fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
%                   fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%         uistack(fl1,'bottom'); 
%         
%         plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);
% 
%         text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    pltHill = plot(hillX, hillY, hillClr,'LineWidth', hillLW);
    pltX0   = plot(0*hillX, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltx1   = plot(0*hillX+9, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltY3H  = plot(hillX, 0*hillY+85/H/1000, hillClr,'LineWidth', hillLW);
    
%     axis([-3 3 -1.4 0]);
    xticks([0 0.5 1 2 3 4 5 6 7 8 9 10]);
%     yticks([-1.4 -1.2 -1.0 -0.8 -0.6 -0.4 -0.2 0]);
%     leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
%     set(leg1,'Visible',legBool);
%     set(leg1,'Color','none');
%     set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'UU_xByd');
end

%% VV_xByD_...

%N = [2,3,3,3,3,3];
N = [2,2,2,2,2,2];

for makePlot = 0 : plotVV_xByd-1

plotOffset = [0.05, 0.5, 1,2,3,4,5,6,7,8]; 
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    yByd(:,line_num) = R_xByd(:, 1 + (line_num-1)*yRCol)/D;
    Mean(:,line_num) = 10*R_xByd(:, 2*xDirIdx+1 + (line_num-1)*yRCol)/Uinf2 + plotOffset(line_num);
%     ub(:,line_num)   = Mean(:,line_num)+N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
%     lb(:,line_num)   = Mean(:,line_num)-N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
%     DET(:,line_num)  = R_xByd_DET(:, xDirIdx + (line_num-1)*(RCol+1))/Uinf2 + plotOffset(line_num);    
end

figure(3)
    hold on;
    grid on;

    ylabel("$y/H$"');
    xlabel("$10v'v' + x/H$");
    
    for line_num=1:line_num_max
        
%         for EXP_num=1:EXP_num_max
%             try
%             plt1(EXP_num) = plot(UU_xByd_EXP{EXP_num,line_num}(:,1), UU_xByd_EXP{EXP_num,line_num}(:,2), ...
%                             EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
%                             'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
%             catch
%             end
%         end
        
        plt2 = plot(Mean(:,line_num),yByd(:,line_num), meanClr,'LineWidth',LW0_5);
%         fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
%                   fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%         uistack(fl1,'bottom'); 
%         
%         plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);
% 
%         text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    pltHill = plot(hillX, hillY, hillClr,'LineWidth', hillLW);
    pltX0   = plot(0*hillX, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltx1   = plot(0*hillX+9, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltY3H  = plot(hillX, 0*hillY+85/H/1000, hillClr,'LineWidth', hillLW);
    
%     axis([-3 3 -1.4 0]);
    xticks([0 0.5 1 2 3 4 5 6 7 8 9 10]);
%     yticks([-1.4 -1.2 -1.0 -0.8 -0.6 -0.4 -0.2 0]);
%     leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
%     set(leg1,'Visible',legBool);
%     set(leg1,'Color','none');
%     set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'VV_xByd');
end

%% UV_xByD_...

%N = [2,3,3,3,3,3];
N = [2,2,2,2,2,2];

for makePlot = 0 : plotUV_xByd-1

plotOffset = [0.05, 0.5, 1,2,3,4,5,6,7,8]; 
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    yByd(:,line_num) = R_xByd(:, 1 + (line_num-1)*yRCol)/D;
    Mean(:,line_num) = -10*R_xByd(:, xDirIdx+1 + (line_num-1)*yRCol)/Uinf2 + plotOffset(line_num);
%     ub(:,line_num)   = Mean(:,line_num)+N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
%     lb(:,line_num)   = Mean(:,line_num)-N(line_num)*R_xByd(:, xDirIdx + RCol + (line_num-1)*yRCol)/Uinf2;
%     DET(:,line_num)  = R_xByd_DET(:, xDirIdx + (line_num-1)*(RCol+1))/Uinf2 + plotOffset(line_num);    
end

figure(4)
    hold on;
    grid on;

    ylabel("$y/H$"');
    xlabel("$-10u'v' + x/H$");
    
    for line_num=1:line_num_max
        
%         for EXP_num=1:EXP_num_max
%             try
%             plt1(EXP_num) = plot(UU_xByd_EXP{EXP_num,line_num}(:,1), UU_xByd_EXP{EXP_num,line_num}(:,2), ...
%                             EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
%                             'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
%             catch
%             end
%         end
        
        plt2 = plot(Mean(:,line_num),yByd(:,line_num), meanClr,'LineWidth',LW0_5);
%         fl1  = fill([yByd', fliplr(yByd')], [ub(:,line_num)', fliplr(lb(:,line_num)')], ...
%                   fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%         uistack(fl1,'bottom'); 
%         
%         plt3 = plot(yByd, DET(:,line_num), DETclr,'LineWidth',LW0_5);
% 
%         text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',txtFSize);
    end
    
    pltHill = plot(hillX, hillY, hillClr,'LineWidth', hillLW);
    pltX0   = plot(0*hillX, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltx1   = plot(0*hillX+9, linspace(1,85/H/1000,Nx), hillClr,'LineWidth', hillLW);
    pltY3H  = plot(hillX, 0*hillY+85/H/1000, hillClr,'LineWidth', hillLW);
    
%     axis([-3 3 -1.4 0]);
    xticks([0 0.5 1 2 3 4 5 6 7 8 9 10]);
%     yticks([-1.4 -1.2 -1.0 -0.8 -0.6 -0.4 -0.2 0]);
%     leg1 = legend([plt1(1) plt1(2) plt1(3) plt3 plt2 fl1],leg_EXP_IPC,'Location','southeast');
%     set(leg1,'Visible',legBool);
%     set(leg1,'Color','none');
%     set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'UV_xByd');
end

%% Save mat file

save([uqPlotsDir '/' caseName '.mat'])


end