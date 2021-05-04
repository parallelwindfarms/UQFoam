function [ ] = periodicHill_main_dataAss(data, caseName, plotU_xByd, plotTauW, plotnut_xByd, tauWTimeDir)
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

sample_xByd = {'0.05','0.50','1.00','2.00','3.00',...
               '4.00','5.00','6.00','7.00','8.00'}; 
sample = {'0_05','0_50','1_00','2_00','3_00',...
          '4_00','5_00','6_00','7_00','8_00'};   
line_num_max  = length(sample);

% DNS data
files = dir([data '0DNS/DNS*']);
DNS = {length(files)};
for i=1:length(files)
    tmpDNS = dlmread([files(i).folder '/' files(i).name], ' ', 37, 0);
    DNS(i) = {tmpDNS};
end
DNStauW = dlmread([data '0DNS/tauW'], ',', 0, 0);

% DET & UQ data
sampleTimeDir = './postProcessing/sample/latestTimeDir';
files = dir([sampleTimeDir '/*xByd']);  
for i=1:length(files)
    load([sampleTimeDir '/' files(i).name]);
end
files = dir([sampleTimeDir '/*c']);
for i=1:length(files)
    load([sampleTimeDir '/' files(i).name]);
end
load('UtauW_dataAss.mat')

% Reading tauW 
sampleTimeDir = ['postProcessing/sampleTauW/surface/' tauWTimeDir];
% files = dir([sampleTimeDir '/tauW*']);
tauWDET   = dlmread([sampleTimeDir '/tauWSpAvg_patch_hills.raw'], ' ', 2, 0);
tauW0     = dlmread([sampleTimeDir '/tauW0SpAvg_patch_hills.raw'], ' ', 2, 0);
tauWSigma = dlmread([sampleTimeDir '/tauWSigmaSpAvg_patch_hills.raw'], ' ', 2, 0);


H    = 28/1000;
D    = H;
R    = D/2;
x    = pc(:,1);
y    = p_xByd(:,1);
tmpxByd = x/D;

% Periodic Hill Geometry
Nx = 500;
[hillX, hillY] = periodicHill_geom(H*1000, Nx);
hillX = hillX/H/1000;
hillY = hillY/H/1000 -1e-3;

hillLS = '-';
hillLW = 1;
hillClr = [.0 .5 .0];

%% Initialization

rhoInf= 1.10;
nu    = 1e-5;
Uinf  = 1;
Uinf2 = Uinf^2;
uCol  = 3;
RCol  = 6;
yuCol = 10;
yRCol = 13;
xDirIdx = 2;
yDirIdx = 3;
yNutCol = 7;

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

DNSclr   = '-k';
DNSedge  = 'none';
DNSsize  = 5;
LW2      = 2;
LW1      = 1;
LW0_5    = 0.5;

DETclr   = [0.6350, 0.0780, 0.1840];%'r';

legBool    = 'on';
leg_DNS_IPC = {'DNS', 'DET', '$\mathbf{E}[\bullet]$', '$\mathbf{E}[\bullet] \pm 2\sqrt{\mathbf{V}[\bullet]}$'};

fSize = 17;
txtFSize = 17;
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', fSize);
% set(0,'DefaultLegendFontSize',fSize);
% set(0,'DefaultTextFontSize',fSize);

% set(groot, 'units', 'inches', 'position', [0 0 8 4])
set(groot, 'defaultFigureUnits','inches')
set(groot, 'defaultFigurePosition',[2.5 1.5 9 3])

% set(groot, 'defaultFigurePaperPositionMode', 'manual');
set(groot, 'defaultFigurePaperUnits', 'inches');
set(groot, 'defaultFigurePaperPosition', [2.5 1.5 9 3]);

%% U_xByD_...

scaling = 2;
N = 2*scaling;

for makePlot = 0 : plotU_xByd-1

plotOffset = [0.05, 0.5, 1,2,3,4,5,6,7,8];  
textOffset = [4.75,3.2,2.2,1.15,0.85,0.28];

for line_num=3:line_num_max
    yByd(:,line_num) = U_xByd(:, 1 + (line_num-1)*yuCol)/D;
    
    Mean(:,line_num) = scaling*U_xByd(:, xDirIdx + (line_num-1)*yuCol)/Uinf + plotOffset(line_num);
    ub(:,line_num)   = Mean(:,line_num)+N*U_xByd(:, xDirIdx+uCol + (line_num-1)*yuCol)/Uinf;
    lb(:,line_num)   = Mean(:,line_num)-N*U_xByd(:, xDirIdx+uCol + (line_num-1)*yuCol)/Uinf;

    Mean_DA(:,line_num) = scaling*U_dataAss(:, line_num)/Uinf + plotOffset(line_num);
    ub_DA(:,line_num)   = Mean_DA(:,line_num)+N*U_Std_dataAss(:, line_num)/Uinf;
    lb_DA(:,line_num)   = Mean_DA(:,line_num)-N*U_Std_dataAss(:, line_num)/Uinf;

end

    
figure(1)
    hold on;
    grid off;

    ylabel("$y/H$"');
    if scaling>1
        xlabel(string(scaling)+"$u/u_{b} + x/H$");
    else
        xlabel("$u/u_{b} + x/H$");
    end
    
    for line_num=3:line_num_max
        
        plt1 = plot(scaling*DNS{line_num}(:,3) + plotOffset(line_num), DNS{line_num}(:,2), ...
                    DNSclr,'MarkerSize',DNSsize,'LineWidth',LW1);
        plt2 = plot(Mean(:,line_num), yByd(:,line_num), 'b','LineWidth',LW1);
        plt3 = plot(Mean_DA(:,line_num), yByH_at_xByH(:,line_num), 'm-','color',[0.6350, 0.0780, 0.1840],'LineWidth',LW1);
        
    end
    
    pltHill = plot(hillX, hillY, hillLS,'LineWidth', hillLW, 'color', hillClr);
    pltX0   = plot(0*hillX, linspace(1,85/H/1000,Nx), hillLS,'LineWidth', hillLW, 'color', hillClr);
    pltx1   = plot(0*hillX+9, linspace(1,85/H/1000,Nx), hillLS,'LineWidth', hillLW, 'color', hillClr);
    pltY3H  = plot(hillX, 0*hillY+85/H/1000, hillLS,'LineWidth', hillLW,'color', hillClr);
    
    axis([0 11 0 3.1]);
    xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
% %     yticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2]);
    leg1 = legend([plt1 plt2 plt3],{'DNS','Prior','Posterior'},'Location','northwest');
    set(leg1,'Visible','off');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'U_xByd_DNS_DataAss',0,1);
    

figure(2)
    hold on;
    grid off;

    ylabel("$y/H$"');
    if scaling>1
        xlabel(string(scaling)+"$u/u_{b} + x/H$");
    else
        xlabel("$u/u_{b} + x/H$");
    end
    
    for line_num=3:line_num_max
        
        plt1 = plot(scaling*DNS{line_num}(:,3) + plotOffset(line_num), DNS{line_num}(:,2), ...
                    DNSclr,'MarkerSize',DNSsize,'LineWidth',LW1);
        fl1  = fill([ub(:,line_num)', fliplr(lb(:,line_num)')], [yByd(:,line_num)', fliplr(yByd(:,line_num)')], ...
                  fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', 'b');
        uistack(fl1,'bottom'); 
        fl2  = fill([ub_DA(:,line_num)', fliplr(lb_DA(:,line_num)')], [yByH_at_xByH(:,line_num)', fliplr(yByH_at_xByH(:,line_num)')], ...
                  [0.6350, 0.0780, 0.1840], 'FaceAlpha', FaceAlpha,'EdgeColor',[0.6350, 0.0780, 0.1840]);
        uistack(fl2,'bottom');         
    end
    
    pltHill = plot(hillX, hillY, hillLS,'LineWidth', hillLW, 'color', hillClr);
    pltX0   = plot(0*hillX, linspace(1,85/H/1000,Nx), hillLS,'LineWidth', hillLW, 'color', hillClr);
    pltx1   = plot(0*hillX+9, linspace(1,85/H/1000,Nx), hillLS,'LineWidth', hillLW, 'color', hillClr);
    pltY3H  = plot(hillX, 0*hillY+85/H/1000, hillLS,'LineWidth', hillLW,'color', hillClr);
    
    axis([0 11 0 3.1]);
    xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
    leg1 = legend([plt1 fl1 fl2],{'DNS','Prior','Posterior'},'Location','northwest');
    set(leg1,'Visible','off');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'U_Std_xByd_DNS_DataAss',0,1);
end

%% tauW_xByD_...

scaling = 1;
N = 2*scaling;

for makePlot = 0 : plotTauW-1

N = 2*scaling;

for makePlot = 0 : plotTauW -1
xByd = tauW0(:,1)/D;
Mean = -scaling*tauW0(:, 4)/Uinf;
ub   = Mean+N*tauWSigma(:, 4)/Uinf;
lb   = Mean-N*tauWSigma(:, 4)/Uinf;

Mean_DA = scaling*tauW_dataAss/Uinf;
ub_DA = Mean_DA+N*tauW_Std_dataAss/Uinf;
lb_DA = Mean_DA-N*tauW_Std_dataAss/Uinf;
    
figure(3)
    hold on;
    grid off;
    
    xlabel("$x/H$");
    if scaling>1
        ylabel(string(scaling)+"$\tau_{w}/u_{b}$");
    else
        ylabel("$\tau_{w}/u_{b}$");
    end

    plt0 = plot(DNStauW(:,1), DNStauW(:,2)*0, 'k--');
    plt1 = plot(DNStauW(:,1), DNStauW(:,2), 'k-','MarkerSize',DNSsize,'LineWidth',LW1);
    plt2 = plot(tauW_xByH, Mean, 'b','LineWidth',LW1);
    plt3 = plot(tauW_xByH, Mean_DA, 'm-','color',[0.6350, 0.0780, 0.1840],'LineWidth',LW1);
        
    axis([0 11 -0.02 0.04]);
    xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
    leg1 = legend([plt1 plt2 plt3],{'DNS', 'Prior', 'Posterior'},'Location','northwest','Orientation','horizontal');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'tauW_xByd_DNS_DataAss',0,2);
    
figure(4)
    hold on;
    grid off;
    
    xlabel("$x/H$");
    if scaling>1
        ylabel(string(scaling)+"$\tau_{w}/u_{b}$");
    else
        ylabel("$\tau_{w}/u_{b}$");
    end

    plt0 = plot(DNStauW(:,1), DNStauW(:,2)*0, 'k--');
    plt1 = plot(DNStauW(:,1), DNStauW(:,2), 'k-','MarkerSize',DNSsize,'LineWidth',LW1);
    fl1  = fill([tauW_xByH', fliplr(tauW_xByH')], [ub', fliplr(lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', 'b');
    uistack(fl1,'bottom'); 
    fl2  = fill([tauW_xByH', fliplr(tauW_xByH')], [ub_DA', fliplr(lb_DA')], ...
              [0.6350, 0.0780, 0.1840], 'FaceAlpha', FaceAlpha,'EdgeColor',[0.6350, 0.0780, 0.1840]);
    uistack(fl2,'bottom');   
        
    axis([0 11 -0.02 0.04]);
    xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
    leg1 = legend([plt1 fl1 fl2],{'DNS', 'Prior', 'Posterior'},'Location','northwest','Orientation','horizontal');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;
    
    uqPlotsPrint(uqPlotsDir,'tauW_Std_xByd_DNS_DataAss',0,2);
    
figure(5)
    hold on;
    grid off;
    
    xlabel("$x_{reattach}/H$");
    
    DNS_xReattch = 5.5;
    DET_xReattch = 3.5;
    
    
    IPC_xReattch = 3.5;
    c = tauWSigma(:, 4)/max(tauWSigma(:, 4));
    c = c(xByd>=2.5 & xByd<=7)';
    IPC_uncBound = linspace(2.5, 7, length(c));
    z = zeros(size(c));
    y = 1.5;
    
    fl1  = surface([IPC_uncBound;IPC_uncBound], [c*0+y;c*0+y], [z;z],...
        [c;c], 'FaceColor','none','EdgeColor','interp',...
        'EdgeAlpha', 0.3,'LineWidth',LW1*10);
    colormap winter;
    plt0 = plot(DNS_xReattch, y, 'sk','MarkerSize',DNSsize,'LineWidth',LW1);
%     plt1 = plot(DET_xReattch, y, '^r','MarkerSize',DNSsize,'LineWidth',LW1);
    plt2 = plot(IPC_xReattch, y, 'ob','MarkerSize',DNSsize,'LineWidth',LW1);
    
    
    DA_xReattch  = 5.8;
    c = tauW_Std_dataAss/max(tauW_Std_dataAss);
    c = c(xByd>=4 & xByd<=6.5)';
    IPC_uncBound = linspace(4, 6.5, length(c));
    z = zeros(size(c));
    y = y-0.75;
    
    fl1  = surface([IPC_uncBound;IPC_uncBound], [c*0+y;c*0+y], [z;z],...
        [c;c], 'FaceColor','none','EdgeColor','interp',...
        'EdgeAlpha', 0.3,'LineWidth',LW1*10);
    colormap winter; 
    plt0 = plot(DNS_xReattch, y, 'sk','MarkerSize',DNSsize,'LineWidth',LW1);
%     plt1 = plot(DET_xReattch, y, '^r','MarkerSize',DNSsize,'LineWidth',LW1);
    plt3 = plot(DA_xReattch, y, 'ob','MarkerSize',DNSsize,'LineWidth',LW1,...
        'color',[0.6350, 0.0780, 0.1840]);
    
    axis([0 11 0 6]);
    xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
    yticks([]);
    leg1 = legend([plt0 plt2 plt3],...
        {'DNS', 'Prior', 'Posterior'},'Location','southwest');
    set(leg1,'Visible',legBool);
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'tauW_xByd_DNS_reattach',0,2);
end

%% Save mat file

% save([uqPlotsDir '/' caseName '.mat'])


end