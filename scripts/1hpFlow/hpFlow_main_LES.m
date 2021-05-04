function [ ] = hpFlow_main_LES(data, caseName, plotp0_c, plotU0_c,...
                plotGradp,plotU_xByd,plotU_xByd_offset,plotU_modes)
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
addpath([SCRIPTS '/1hpFlow']);
addpath([SCRIPTS '/matlabScripts']);
addpath([SCRIPTS '/matlabScripts/matlab2tikz-master/src']);
addpath([SCRIPTS '/matlabScripts/export_fig/export_fig-master']);

load(['./gradP1_P0.mat']);
load(['./gradP2_P0.mat']);
load(['./U_modes_IPC_xByd_2.mat']);
load(['./U_modes_NIPC_xByd_2.mat']);

addpath("./postProcessing/sample/10");

files = [   
            "centerline_pMean0_pMean1_pMeanSigma.xy";
            "centerline_UMean0_UMean1_UMeanSigma.xy";
            "X_by_delta_2_pMean0_pMean1_pMeanSigma.xy";
            "X_by_delta_2_UMean0_UMean1_UMeanSigma.xy";
            "X_by_delta_30_pMean0_pMean1_pMeanSigma.xy";
            "X_by_delta_30_UMean0_UMean1_UMeanSigma.xy";
            "X_by_delta_5_pMean0_pMean1_pMeanSigma.xy";
            "X_by_delta_5_UMean0_UMean1_UMeanSigma.xy"
        ];
    
for idx = 1:length(files)
    load(files(idx));
end

delta       = 0.1;
x           = centerline_pMean0_pMean1_pMeanSigma(:,1);
y           = X_by_delta_2_pMean0_pMean1_pMeanSigma(:,1);
xByd        = x/delta;
yByd        = y/delta;

p0_c        = centerline_pMean0_pMean1_pMeanSigma(:,2);
p1_c        = centerline_pMean0_pMean1_pMeanSigma(:,3);
U0_c        = centerline_UMean0_UMean1_UMeanSigma(:,2);
pSigma_c    = centerline_pMean0_pMean1_pMeanSigma(:,4);
USigma_c    = centerline_UMean0_UMean1_UMeanSigma(:,8);

p0_xByd     = [ X_by_delta_2_pMean0_pMean1_pMeanSigma(:,2) ...
                X_by_delta_5_pMean0_pMean1_pMeanSigma(:,2) ...
                X_by_delta_30_pMean0_pMean1_pMeanSigma(:,2) ];

pSigma_xByd = [ X_by_delta_2_pMean0_pMean1_pMeanSigma(:,4) ...
                X_by_delta_5_pMean0_pMean1_pMeanSigma(:,4) ...
                X_by_delta_30_pMean0_pMean1_pMeanSigma(:,4) ];
            
U0_xByd     = [ X_by_delta_2_UMean0_UMean1_UMeanSigma(:,2) ...
                X_by_delta_5_UMean0_UMean1_UMeanSigma(:,2) ...
                X_by_delta_30_UMean0_UMean1_UMeanSigma(:,2) ];
        
USigma_xByd = [ X_by_delta_2_UMean0_UMean1_UMeanSigma(:,8) ...
                X_by_delta_5_UMean0_UMean1_UMeanSigma(:,8) ...
                X_by_delta_30_UMean0_UMean1_UMeanSigma(:,8) ];
            

%% Initialization

Ub		= 1;
nu0		= 2e-3;
nu1		= 4e-4;
nu2		= 16e-8;

%% Plot settings

% Defaults 
uqClr_NIPC  = {rgb('red'),rgb('red'),rgb('red'),rgb('red'),rgb('red')};
uqClr_IPC = {rgb('blue'),rgb('blue'),rgb('blue'),rgb('blue'),rgb('blue')};
uqLS   = {'-','--','--','-.',':'};
uqMark = {'none','none','none','none','none'};

anaClr   = '-k';
LW1     = 0.25;
LW0_5   = 0.25;

meanClr = '-b';
fill_color = rgb('blue');%[.5 .5 .5];
FaceAlpha = 0.2;
EdgeColor = 'none';

fSize = 25;
txtFSize = 20;
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', fSize);
% set(0,'DefaultLegendFontSize',fSize);
% set(0,'DefaultTextFontSize',fSize);

% set(groot, 'units', 'inches', 'position', [0 0 8 4])
set(groot, 'defaultFigureUnits','inches');
% set(groot, 'defaultFigurePosition',[2.5 1.5 8 5]); % Poster
set(groot, 'defaultFigurePosition',[2.5 1.5 8 6]);

% set(groot, 'defaultFigurePaperPositionMode', 'manual');
set(groot, 'defaultFigurePaperUnits', 'inches');
% set(groot, 'defaultFigurePaperPosition', [2.5 1.5 8 5]); % Poster
set(groot, 'defaultFigurePaperPosition', [2.5 1.5 8 6]);

%% p0_c
N = 2                    % times sigma

for makePlot = 0 : plotp0_c-1
slope = -0.7212/(5-3.798);
p_in  = 0.7212-slope*3.798;

mean = p0_c/p_in;
ub = mean+N*pSigma_c/p_in;
lb = mean-N*pSigma_c/p_in;

figure(1)
    hold on;
    grid off;
    xlabel("$x/\delta$"');
    ylabel("$p/p_{in}$");
   
    plt1 = plot(xByd, (p_in+slope*xByd*delta)/p_in, anaClr,'LineWidth',LW1);
    
    plt2 = plot(xByd, mean, meanClr,'LineWidth',LW0_5);
    fl1 = fill([xByd', fliplr(xByd')], [ub', fliplr(lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%     uistack(fl1,'bottom');

    axis([0 50 0 2])
    xticks([0 10 20 30 40 50]);
    leg1 = legend([plt1 plt2 fl1],{'$1-\frac{x}{p_{in}} \frac{dp}{dx}$', ...
                           '$\mathbf{E}[p]$','$\pm 2\sqrt{\mathbf{V}[p]}$'},'location','northeast');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'p0_c');
end

%% U0_c
for makePlot = 0 : plotU0_c-1
u_ana  = 1.5*ones(length(xByd),1)/Ub;

mean = U0_c/Ub;
ub = mean+N*USigma_c/Ub;
lb = mean-N*USigma_c/Ub;

figure(2)
    hold on;
    grid off;

    xlabel("$x/\delta$"');
    ylabel("$u/u_{avg}$");
    
    plt1 = plot(xByd, u_ana, anaClr,'LineWidth',LW1);
    
    plt2 = plot(xByd, mean, meanClr,'LineWidth',LW0_5);
    fl1 = fill([xByd', fliplr(xByd')], [ub', fliplr(lb')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%     uistack(fl1,'bottom');
    axis([0 50 1 1.6])
%     yticks([1 1.1 1.2 1.3 1.4 1.5]);
    xticks([0 10 20 30 40 50]);

    leg1 = legend([plt1 plt2 fl1],{'$u_{max}/u_{avg}$','$\mathbf{E}[u]$','$\pm 2\sqrt{\mathbf{V}[u]}$'},'location','southeast');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'U0_c');
end

%% p_grad
for makePlot = 0 : plotGradp-1
localFac1 = 10;
localFac2 = 100;
xBydVal = linspace(0,50,length(gradP1_P0));
anaVal  = [(nu1/nu0)*ones(length(gradP1_P0),1)*localFac1 (nu2/nu0)*ones(length(gradP2_P0),1)*localFac2];

mean = [gradP1_P0*localFac1 0.4*gradP2_P0*localFac2];

figure(3)
    hold on;
    grid off;

    xlabel("$x/\delta$"');
    ylabel("$(\partial p_i / \partial x)/(\partial p_0 / \partial x)$");
    
    plt1 = plot(xBydVal, anaVal(:,1), anaClr,'LineWidth',1); 
    plt2 = plot(xBydVal, mean(:,1), meanClr,'LineWidth',1);

    plt3 = plot(xBydVal, anaVal(:,2), '--k','LineWidth',1); 
    plt4 = plot(xBydVal, mean(:,2), '--b','LineWidth',1);    
        
    axis([0 50 -0.5 2]);
    yticks(2*[-0.5 -0.25 0.0 0.25 0.5 0.75 1]);
    xticks([0 10 20 30 40 50]);
        
    leg1 = legend([plt1 plt3 plt2 plt4], ...
        {'$\ \nu_1/ \nu_0 \ \times 10^1$','$\ \nu_2 / \nu_0 \ \times 10^2$',...
        '$\frac{\partial p_1}{\partial x}/\frac{\partial p_0}{\partial x} \times 10^1$',...
        '$\frac{\partial p_2}{\partial x}/\frac{\partial p_0}{\partial x} \times 10^2$'},'location','northeast');
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,'gradP_c');
end

%% U_xByd_2,5,30
for makePlot = 0 : plotU_xByd-1
u_ana  = 1.5*(1-4*(yByd/2).^2)/Ub;

mean = U0_xByd/Ub;
ub = mean+N*USigma_xByd/Ub;
lb = mean-N*USigma_xByd/Ub;

for i=1:3
figure(3 + i)
    hold on;
    grid off;

    xlabel("$y/\delta$"');
    ylabel("$u/u_{avg}$");
    
    plt1 = plot(yByd, u_ana, anaClr,'LineWidth',LW1);
    
    xBydi = i;
    plt2 = plot(yByd, mean(:,xBydi), meanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], [ub(:,xBydi)', fliplr(lb(:,xBydi)')], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
    uistack(fl1,'bottom');
    axis([0 1 0 1.6])
    xticks([0 0.2 0.4 0.6 0.8 1.0]);

    leg1 = legend([plt1 plt2 fl1],{'[$u_{HP}]_{\nu=\mathbf{E}[\nu]}$','$\mathbf{E}[u]$','$\pm 2\sqrt{\mathbf{V}[u]}$'});
    set(leg1,'Color','none');
    set(leg1,'EdgeColor','none');
    hold off;

    uqPlotsPrint(uqPlotsDir,['U0_xByd' num2str(i)]);
end
    
end

%% U_xByd_2_5_30 with offset
for makePlot = 0 : plotU_xByd_offset-1
u_ana  = 1.5*(1-4*(yByd/2).^2)/Ub;

mean = U0_xByd/Ub;
ub = mean+N*USigma_xByd/Ub;
lb = mean-N*USigma_xByd/Ub;

figure(7)
    hold on;
    grid off;
    
    xlabel("$y/\delta$"');
    ylabel("$u/u_{avg}$");
    
    offset = 0.5;
   
for i=1:3
    plt1 = plot(yByd, u_ana + (i-1)*offset, anaClr,'LineWidth',LW1);
    
    xBydi = i;
    plt2 = plot(yByd, mean(:,xBydi) + (i-1)*offset, meanClr,'LineWidth',LW0_5);
    fl1 = fill([yByd', fliplr(yByd')], ...
               [ub(:,xBydi)'+ (i-1)*offset, fliplr(lb(:,xBydi)'+ (i-1)*offset)], ...
              fill_color, 'FaceAlpha', FaceAlpha,'EdgeColor', EdgeColor);
%     uistack(fl1,'bottom');

    if i==1
        text(0.025, 1.5 + (i-1)*offset + 0.1,'$x/\delta = 2$','FontSize',txtFSize)
    elseif i==2
        text(0.025, 1.5 + (i-1)*offset + 0.1,'$x/\delta = 5$','FontSize',txtFSize)
    else
        text(0.025, 1.5 + (i-1)*offset + 0.1,'$x/\delta = 20$','FontSize',txtFSize)
    end        
    
end
    axis([0 1 0 3])  
    yticks([0 0.5 1.0 1.5 2.0 2.5 3.0]);
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    leg1 = legend([plt1 plt2 fl1],{'[$u_{HP}]_{\nu=\mathbf{E}[\nu]}$','$\mathbf{E}[u]$','$\pm 2\sqrt{\mathbf{V}[u]}$'});
    set(leg1,'Color','none','EdgeColor','none');
    hold off;
    uqPlotsPrint(uqPlotsDir,'U0_xByd_offset');
end

%% U_modes

for makePlot = 0 : plotU_modes-1

figure(8)
    hold on;
    grid off;
    
    xlabel("$y/\delta$"');
    ylabel("$u_i/u_{avg}$");
   
    for k = 2 : 3   
    plot(U_modes_IPC_xByd_2(:,1)/delta, U_modes_IPC_xByd_2(:,k+1) * 100, ...
        'color', 'k','LineWidth', 1,'LineStyle',uqLS{k-1},'marker',uqMark{k-1})
    end  
    for k = 2 : 3   
        plot(U_modes_IPC_xByd_2(:,1)/delta, U_modes_IPC_xByd_2(:,k+1) * 100, ...
            'color', uqClr_IPC{k-1},'LineWidth', 1,'LineStyle',uqLS{k-1},'marker',uqMark{k-1})
        plot(U_modes_NIPC_xByd_2(:,1)/delta, U_modes_NIPC_xByd_2(:,k+1) * 100, ...
            'color', uqClr_NIPC{k-1},'LineWidth', 1,'LineStyle',uqLS{k-1},'marker',uqMark{k-1})
    end
   
    axis([0 1 -4 4])  
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    leg1 = legend('$u_1 \times 10^2$','$u_2 \times 10^2$');
    set(leg1,'Color','none','EdgeColor','none');
    hold off;
    uqPlotsPrint(uqPlotsDir,'U_modes_12');

figure(9)   
    hold on;
    grid off;
    
    xlabel("$y/\delta$"');
    ylabel("$u_i/u_{avg}$");
   
    for k = 4 : (size(U_modes_IPC_xByd_2,2)-1)   
    plot(U_modes_IPC_xByd_2(:,1)/delta, U_modes_IPC_xByd_2(:,k+1) * 1000, ...
        'color', 'k','LineWidth', 1,'LineStyle',uqLS{k-1},'marker',uqMark{k-1})
    end  
    for k = 4 : (size(U_modes_IPC_xByd_2,2)-1)  
        plot(U_modes_IPC_xByd_2(:,1)/delta, U_modes_IPC_xByd_2(:,k+1) * 1000, ...
            'color', uqClr_IPC{k-1},'LineWidth', 1,'LineStyle',uqLS{k-1},'marker',uqMark{k-1})
        plot(U_modes_NIPC_xByd_2(:,1)/delta, U_modes_NIPC_xByd_2(:,k+1) * 1000, ...
            'color', uqClr_NIPC{k-1},'LineWidth', 1,'LineStyle',uqLS{k-1},'marker',uqMark{k-1})
    end
   
    axis([0 1 -3 3]) 
    xticks([0 0.2 0.4 0.6 0.8 1.0]);
    leg1 = legend('$u_3 \times 10^3$','$u_4 \times 10^3$','$u_5 \times 10^3$');
    set(leg1,'Color','none','EdgeColor','none');
    hold off;
    uqPlotsPrint(uqPlotsDir,'U_modes_345');
    
%% Save mat file
save([uqPlotsDir '/' caseName '.mat'])

end