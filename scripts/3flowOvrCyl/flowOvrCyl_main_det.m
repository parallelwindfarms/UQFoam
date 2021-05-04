function [ ] = flowOvrCyl_main_det(data, caseName, plotU_c, plotU_xByd, plotV_xByd, plotUU_xByd, plotVV_xByd, plotCd, plotCl, plotCp)
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

% load([uqPlotsDir '/forceCoeffs.mat']);
% t        = coefficient200(:,1);
% Cd       = coefficient200(:,3);
% Cd_mean  = coefficient200_mean(:,3);
% Cl       = coefficient200(:,4);
% Cl_mean  = coefficient200_mean(:,4);

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
yuCol = 4;
yRCol = 7;
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

legBool    = 'off';
leg_EXP = {'EXP 1','EXP 2','EXP 3','DET'};

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

%% U_c

for makePlot = 0 : plotU_c-1

Mean = Uc(:,2)/Uinf;

figure(1)
    hold on;
    grid on;

    xlabel("$x/D$"');
    ylabel("$u / U_{\infty}$");
    
    for EXP_num=1:EXP_num_max
        plt1(EXP_num) = plot(Uc_EXP{EXP_num}(:,1), Uc_EXP{EXP_num}(:,2), ...
                        EXPclr{EXP_num},'MarkerFaceColor',EXPfill{EXP_num}, ...
                        'MarkerSize',EXPsize,'MarkerEdgeColor',EXPedge{EXP_num});
    end 
    
    plt2 = plot(xByd, Mean, meanClr,'LineWidth',LW0_5);
    axis([0 10 -0.4 1])

    leg1 = legend([plt1(1) plt1(2) plt1(3) plt2],leg_EXP,'Location','southeast');
    set(leg1,'Visible','on');
    hold off;

    uqPlotsPrint(uqPlotsDir,'U_centerline','0','1');
end

clear Mean; clear ub; clear lb;

%% U_xByD_1.06,1.54,2.02,4,7,10

for makePlot = 0 : plotU_xByd-1

plotOffset = [4.5,3.0,2.0,1.0,0.75,0.2]-1;  
textOffset = [4.75,3.2,2.2,1.15,0.85,0.28];

for line_num=1:line_num_max
    Mean(:,line_num) = U_xByd(:, xDirIdx + (line_num-1)*yuCol)/Uinf + plotOffset(line_num);
end

figure(2)
    hold on;
    grid on;

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
        
        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',8);
    end
    
    axis([-3 3 -1 5])
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt2],leg_EXP,'Location','southeast');
    set(leg1,'Visible',legBool);    
    hold off;

    uqPlotsPrint(uqPlotsDir,'U_xByd');
end

%% V_xByD_1.06,1.54,2.02,4,7,10

for makePlot = 0 : plotV_xByd-1

plotOffset = [-0.7,-1.4,-2.1,-2.5,-2.75,-2.9];  
textOffset = plotOffset + 0.075;

for line_num=1:line_num_max
    Mean(:,line_num) = U_xByd(:, yDirIdx + (line_num-1)*yuCol)/Uinf + plotOffset(line_num);
end

figure(3)
    hold on;
    grid on;

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
        
        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',8);
    end
    
    axis([-3 3 -3 0])
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt2 ],leg_EXP,'Location','southeast');
    set(leg1,'Visible',legBool);    
    hold off;

    uqPlotsPrint(uqPlotsDir,'V_xByd');
end

%% UU_xByD_1.06,1.54,2.02,4,7,10
 
for makePlot = 0 : plotUU_xByd-1

plotOffset = [-0.3,-0.55,-0.8,-0.95,-1.1,-1.25];  
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    Mean(:,line_num) = R_xByd(:, xDirIdx + (line_num-1)*yRCol)/Uinf2 + plotOffset(line_num);
end

figure(4)
    hold on;
    grid on;

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
        
        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',8);
    end
    
    xlim([-3 3]);
    ylim([-1.4 0]);
%     axis([-3 3 -1.5 0])
%     yticks([-1.6 -1.2 -0.8 -0.4 0]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt2],leg_EXP,'Location','southeast');
    set(leg1,'Visible',legBool);    
    hold off;

    uqPlotsPrint(uqPlotsDir,'UU_xByd');
end

%% VV_xByD_1.06,1.54,2.02,4,7,10

factor = [1.5,1.5,1,1,1,1];

for makePlot = 0 : plotVV_xByd-1

EXPOffset = [0.2,0,0,0,0,0];
plotOffset = [-0.4,-0.9,-1.4,-1.6,-1.75,-2.05];  
textOffset = plotOffset + 0.05;

for line_num=1:line_num_max
    Mean(:,line_num) = R_xByd(:, 2*xDirIdx+1 + (line_num-1)*yRCol)/Uinf2/factor(line_num) + plotOffset(line_num);
end

figure(5)
    hold on;
    grid on;

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
        
        text(2, textOffset(line_num),['$x/D =$ ' sample_xByd{line_num}],'FontSize',8);
    end
    
    xlim([-3 3]);
%     axis([-3 3 -2.25 0])
%     yticks([-2.25 -1.75 -1.25 -0.75 -0.25]);
    leg1 = legend([plt1(1) plt1(2) plt1(3) plt2],leg_EXP,'Location','southeast');
    set(leg1,'Visible',legBool);    
    hold off;

    uqPlotsPrint(uqPlotsDir,'VV_xByd');
end

%% Cd

for makePlot = 0 : plotCd-1

TUinf_D = t*Uinf/D + 600;

figure(6)
    hold on;
    grid on;

    xlabel("$tU_{\infty}/D$"');
    ylabel("$C_d$");
      
    Mean = Cd;
    plt1 = plot(TUinf_D, Mean, tclr,'LineWidth',LW1);
        
    Mean = Cd_mean;
    plt2 = plot(TUinf_D, Mean, meanClr,'LineWidth',LW1);

    ylim([1 1.6]);
%     yticks([1 1.1 1.2 1.3 1.4 1.5]);
    xlim([1200 1700]);
    xticks(linspace(1200,1700,6));

    leg1 = legend([plt1 plt2],{'$C_d$','$<C_d>$'},'Location','northeast');
    set(leg1,'Visible','on');
    hold off;

    fprintf('Cd_mean = %f\n', mean(Cd_mean(:,1)))
    
    set(figure(6), 'Units','inches')
    set(figure(6), 'Position',[2.5 1.5 4 3])
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [2.5 1.5 4 3]);
    uqPlotsPrint(uqPlotsDir,'Cd');
end

%% Cl

for makePlot = 0 : plotCl-1

TUinf_D = t*Uinf/D + 600;

figure(7)
    hold on;
    grid on;

    xlabel("$tU_{\infty}/D$"');
    ylabel("$C_l$");
      
    Mean = Cl;
    plt1 = plot(TUinf_D, Mean, tclr,'LineWidth',LW1);

    
    Mean = Cl_mean;
    plt2 = plot(TUinf_D, Mean, meanClr,'LineWidth',LW1);
    
    ylim([-1 3]);
    xlim([1200 1700]);
    xticks(linspace(1200,1700,6));

    leg1 = legend([plt1 plt2],{'$C_l$','$<C_l>$'},'Location','northeast');
    hold off;

    fprintf('Cl_mean = %f\n', mean(Cl_mean(:,1)))
    
    set(figure(7), 'Units','inches')
    set(figure(7), 'Position',[2.5 1.5 4 3])
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [2.5 1.5 4 3]);
    uqPlotsPrint(uqPlotsDir,'Cl');
end

%% Cp

for makePlot = 0 : plotCp-1
  
p_mean = load([uqPlotsDir '/p_mean_cyl.csv']);

Cp_factor = 2/(rhoInf*Uinf2);
Mean  = Cp_factor*p_mean(:,2);
theta = 180*(1-p_mean(:,1)/max(p_mean(:,1)));

figure(8)
    hold on;
    grid on;

    xlabel("$\theta$"');
    ylabel("$C_p$");

    plt1 = plot(Cp_EXP(:,1), Cp_EXP(:,2), 'ks','MarkerFaceColor',...
        EXPfill{1}, 'MarkerSize',7.5,'MarkerEdgeColor',EXPedge{1});
    
    plt2 = plot(theta, Mean, meanClr,'LineWidth',LW1);

    ylim([-1.5 1.5]);
    yticks([-1.5 -1 -0.5 0 0.5 1 1.5]);
    xlim([0 180]);
    xticks([0 30 60 90 120 150 180]);

    leg1 = legend([plt1 plt2],{'EXP','DET'}, 'Location','northeast');
    set(leg1,'Visible','on');
    hold off;
    
    set(figure(8), 'Units','inches')
    set(figure(8), 'Position',[2.5 1.5 8 3])
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [2.5 1.5 8 3]);
    uqPlotsPrint(uqPlotsDir,'Cp');
end

%% Save mat file

save([uqPlotsDir '/' caseName '.mat'])


end