%% Clean-up

clear all;
clc;
close all;

if exist("uqPlots", 'dir')
    rmdir("uqPlots", 's');
end
mkdir uqPlots;

plotp    = 0;
plotU    = 0;
plotUrms = 0;
plotuv   = 0;

plotp    = 1;
plotU    = 1;
plotUrms = 1;
plotuv   = 1;

%% Load data

addpath("/home/p285464/scripts/DNS");               % cluster data
addpath("/home/jigar/Peregrine/data/scripts/DNS");  % mounted data

Udns = dlmread("UMean.txt", '', 2, 0);
Rdns = dlmread("RMean.txt", '', 2, 0);

addpath("./postProcessing/collapsedFields/latesTimeDir");

files = [   
    
        "p0Mean.xy"; "pSigmaMean.xy";
        "pMean0.xy"; "pMeanSigma.xy";
        
        "U0Mean_X.xy"; "USigmaMean_X.xy"; 
        "U0Prime2Mean_XX.xy"; "U0Prime2Mean_YY.xy"; "U0Prime2Mean_ZZ.xy";"U0Prime2Mean_XY.xy";
        "USigmaPrime2Mean_XX.xy"; "USigmaPrime2Mean_YY.xy"; "USigmaPrime2Mean_ZZ.xy";"USigmaPrime2Mean_XY.xy";
        
        "UMean0_X.xy"; "UMeanSigma_X.xy";
        "RMean0_XX.xy"; "RMean0_YY.xy"; "RMean0_ZZ.xy"; "RMean0_XY.xy";
        "RMeanSigma_XX.xy"; "RMeanSigma_YY.xy"; "RMeanSigma_ZZ.xy"; "RMeanSigma_XY.xy";
        
        "UMeanSigma2_XX.xy"; "UMeanSigma2_YY.xy"; "UMeanSigma2_ZZ.xy"; "UMeanSigma2_XY.xy";
        "UMeanSigmaPrime2Mean_XX.xy"; "UMeanSigmaPrime2Mean_YY.xy"; "UMeanSigmaPrime2Mean_ZZ.xy"; "UMeanSigmaPrime2Mean_XY.xy";
        
        ];
    
for idx = 1:length(files)
    load(files(idx));
end

%% Initialization

col1 = 1;
col2 = 2;

M     = dlmread("postProcessing/patchExpression_yPlus/0/bottomWall", '', 1,0);
YPlus = mean(M(:,col2))
y     = p0Mean(:,1);
yWall = y(2);

Ub		= 0.1335;
nu		= 2e-5;


UT		= 0.0079;               % DNS @ Re_t = 395
UT2		= UT*UT;
UT_nu	= UT/nu;
Y       = Udns(:,1);
Ypl     = Udns(:,2);
N		= 2;                    % times sigma

%% pMean

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

p = pMean0(:,col2);
pSigmaMean = pSigmaMean(:,col2);
pMeanSigma = pMeanSigma(:,col2);

if (plotp)

figure(11)
hold on;
grid on;
xlabel("y/\delta");
ylabel("<p>")
axis([0 1 -5e-4 5e-4])
plot(y, p, 'r');
plot(y, p + N*pSigmaMean, 'b', y, p - N*pSigmaMean, 'b');
plot(y, p + N*pMeanSigma, 'g', y, p - N*pMeanSigma, 'g');
% legend('U','\pm sigmaMean','','\pm meanSigma','')
hold off;
saveas(figure(11),'./uqPlots/0pMean_vs_y', 'epsc')

figure(22)
semilogx(yPlus, p, 'r');
hold on;
grid on;
xlabel("y^+");
ylabel("<p>")
axis([0 400 -5e-4 5e-4])
semilogx(yPlus, p + N*pSigmaMean, 'b', yPlus, p - N*pSigmaMean, 'b');
semilogx(yPlus, p + N*pMeanSigma, 'g', yPlus, p - N*pMeanSigma, 'g');
% legend('U','\pm sigmaMean','','\pm meanSigma','')
hold off;
saveas(figure(22),'./uqPlots/0pMean_vs_yPl', 'epsc')

end

%% UMean

ut      = 1.8*nu/yWall 
% ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

U = U0Mean_X(:,col2)/ut;
USigmaMean = USigmaMean_X(:,col2)/ut;
UMeanSigma = UMeanSigma_X(:,col2)/ut;

if (plotU)

figure(1)
plot(Y, Udns(:,3), 'k');
hold on;
grid on;
xlabel("y/\delta");
ylabel("<u>/u_\tau")
axis([0 1 0 30])
plot(y, U, 'r');
plot(y, U + N*USigmaMean, 'b', y, U - N*USigmaMean, 'b');
plot(y, U + N*UMeanSigma, 'g', y, U - N*UMeanSigma, 'g');
% legend('U','\pm sigmaMean','','\pm meanSigma','')
hold off;
saveas(figure(1),'./uqPlots/1UMean_vs_y', 'epsc')

figure(2)
semilogx(Ypl, Udns(:,3), 'k');
hold on;
grid on;
xlabel("y^+");
ylabel("<u>/u_\tau")
axis([0 400 0 30])
semilogx(yPlus, U, 'r');
semilogx(yPlus, U + N*USigmaMean, 'b', yPlus, U - N*USigmaMean, 'b');
semilogx(yPlus, U + N*UMeanSigma, 'g', yPlus, U - N*UMeanSigma, 'g');
% legend('U','\pm sigmaMean','','\pm meanSigma','')
hold off;
saveas(figure(2),'./uqPlots/1UMean_vs_yPl', 'epsc')

end

%% UU

% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

U0Prime2Mean = [(U0Prime2Mean_XX(:,col2)), (U0Prime2Mean_YY(:,col2)), (U0Prime2Mean_ZZ(:,col2))]/ut2;
USigmaPrime2Mean = [(USigmaPrime2Mean_XX(:,col2)), (USigmaPrime2Mean_YY(:,col2)), (USigmaPrime2Mean_ZZ(:,col2))]/ut2;

RMean0 = [(RMean0_XX(:,col2)), (RMean0_YY(:,col2)), (RMean0_ZZ(:,col2))]/ut2;
RMeanSigma = [(RMeanSigma_XX(:,col2)), (RMeanSigma_YY(:,col2)), (RMeanSigma_ZZ(:,col2))]/ut2;

UMeanSigma2 = [(UMeanSigma2_XX(:,col2)), (UMeanSigma2_YY(:,col2)), (UMeanSigma2_ZZ(:,col2))]/ut2;
UMeanSigmaPrime2Mean = [(UMeanSigmaPrime2Mean_XX(:,col2)), (UMeanSigmaPrime2Mean_YY(:,col2)), (UMeanSigmaPrime2Mean_ZZ(:,col2))]/ut2;

if (plotUrms)
    
offset = 1;
ymax = 6*offset;

figure(3)
hold on;
grid on;
xlabel("y/\delta");
ylabel("( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau")
axis([0 1 0 ymax])
tmp = sqrt(Rdns(:,3:5));
plot(Y, tmp(:,2), 'k', Y, tmp(:,3) + 1*offset, 'k', Y, tmp(:,1) + 2*offset, 'k');
tmp = sqrt(U0Prime2Mean);
plot(y, tmp(:,2), 'r', y, tmp(:,3) + 1*offset, 'r', y, tmp(:,1) + 2*offset, 'r');
tmp = sqrt(U0Prime2Mean + N*USigmaPrime2Mean);
plot(y, tmp(:,2), 'g', y, tmp(:,3) + 1*offset, 'g', y, tmp(:,1) + 2*offset, 'g');
tmp = sqrt(U0Prime2Mean - N*USigmaPrime2Mean);
plot(y, tmp(:,2), 'g', y, tmp(:,3) + 1*offset, 'g', y, tmp(:,1) + 2*offset, 'g');
hold off;
saveas(figure(3),'./uqPlots/2USigmaPrime2Mean_vs_y', 'epsc')

figure(4)
hold on;
grid on;
xlabel("y/\delta");
ylabel("( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau")
axis([0 1 0 ymax])
tmp = sqrt(Rdns(:,3:5));
plot(Y, tmp(:,2), 'k', Y, tmp(:,3) + 1*offset, 'k', Y, tmp(:,1) + 2*offset, 'k');
tmp = sqrt(RMean0);
plot(y, tmp(:,2), 'r', y, tmp(:,3) + 1*offset, 'r', y, tmp(:,1) + 2*offset, 'r');
tmp = sqrt(RMean0 + N*RMeanSigma);
plot(y, tmp(:,2), 'g', y, tmp(:,3) + 1*offset, 'g', y, tmp(:,1) + 2*offset, 'g');
tmp = sqrt((RMean0 - N*RMeanSigma));
plot(y, tmp(:,2), 'g', y, tmp(:,3) + 1*offset, 'g', y, tmp(:,1) + 2*offset, 'g');
hold off;
saveas(figure(4),'./uqPlots/3RMeanSigma_vs_y', 'epsc')

figure(5)
hold on;
grid on;
xlabel("y/\delta");
ylabel("( v_{rms} , w_{rms} + u_\tau , u_{rms} + 2u_\tau ) / u_\tau")
axis([0 1 0 ymax])
tmp = sqrt(Rdns(:,3:5));
plot(Y, tmp(:,2), 'k', Y, tmp(:,3) + 1*offset, 'k', Y, tmp(:,1) + 2*offset, 'k');
tmp = sqrt(RMean0);
plot(y, tmp(:,2), 'r', y, tmp(:,3) + 1*offset, 'r', y, tmp(:,1) + 2*offset, 'r');
tmp = sqrt(RMean0 + N*USigmaPrime2Mean);
plot(y, tmp(:,2), 'g', y, tmp(:,3) + 1*offset, 'g', y, tmp(:,1) + 2*offset, 'g');
tmp = sqrt(RMean0 - N*USigmaPrime2Mean);
plot(y, tmp(:,2), 'g', y, tmp(:,3) + 1*offset, 'g', y, tmp(:,1) + 2*offset, 'g');
hold off;
saveas(figure(5),'./uqPlots/4R0USigmaPrime2Mean_vs_y', 'epsc')

end

%% uv
    
% ut      = 1.8*nu/yWall 
ut		= YPlus*nu/yWall 		% LES
ut2		= ut*ut;
ut_nu 	= ut/nu;
yPlus   = y*ut_nu;

U0Prime2Mean = [U0Prime2Mean, -U0Prime2Mean_XY(:,col2)/ut2];
USigmaPrime2Mean = [USigmaPrime2Mean, USigmaPrime2Mean_XY(:,col2)/ut2];

RMean0 = [RMean0, -RMean0_XY(:,col2)/ut2];
RMeanSigma = [RMeanSigma, RMeanSigma_XY(:,col2)/ut2];

UMeanSigma2 = [UMeanSigma2, UMeanSigma2_XY(:,col2)/ut2];
UMeanSigmaPrime2Mean = [UMeanSigmaPrime2Mean, UMeanSigmaPrime2Mean_XY(:,col2)/ut2];

if(plotuv)

figure(7)
hold on;
grid on;
xlabel("y/\delta");
ylabel("-uv/u_\tau^2")
axis([0 1 0 1])
tmp = -Rdns(:,6);
plot(Y, tmp(:,1), 'k');
tmp = U0Prime2Mean;
plot(y, tmp(:,4), 'r');
tmp = U0Prime2Mean + N*USigmaPrime2Mean;
plot(y, tmp(:,4), 'g');
tmp = U0Prime2Mean - N*USigmaPrime2Mean;
plot(y, tmp(:,4), 'g');
hold off;
saveas(figure(7),'./uqPlots/2uvUSigmaPrime2Mean_vs_y', 'epsc')

figure(8)
hold on;
grid on;
xlabel("y/\delta");
ylabel("-uv/u_\tau^2")
axis([0 1 0 1])
tmp = -Rdns(:,6);
plot(Y, tmp(:,1), 'k');
tmp = RMean0;
plot(y, tmp(:,4), 'r');
tmp = RMean0 + N*RMeanSigma;
plot(y, tmp(:,4), 'g');
tmp = RMean0 - N*RMeanSigma;
plot(y, tmp(:,4), 'g');
hold off;
saveas(figure(8),'./uqPlots/3uvRMeanSigma_vs_y', 'epsc')

figure(9)
hold on;
grid on;
xlabel("y/\delta");
ylabel("-uv/u_\tau^2")
axis([0 1 0 1])
tmp = -Rdns(:,6);
plot(Y, tmp(:,1), 'k');
tmp = RMean0;
plot(y, tmp(:,4), 'r');
tmp = RMean0 + N*USigmaPrime2Mean;
plot(y, tmp(:,4), 'g');
tmp = RMean0 - N*USigmaPrime2Mean;
plot(y, tmp(:,4), 'g');
hold off;
saveas(figure(9),'./uqPlots/4uvR0USigmaPrime2Mean_vs_y', 'epsc')


end