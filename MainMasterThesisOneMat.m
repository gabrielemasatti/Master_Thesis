% Msc Thesis Masatti Gabriele
%
% Based on Patrizia Semeraro:
% Multivariate tempered stable additive subordination for financial models
%
% Supervisors: Michele Azzone, Roberto Baviera
% 
% -*- -*- -*- -*- -*- -*- -*- -*- -*-
%
% Academic Year 2022 - 2023
%
% RunMain

%% clearing
clear
close all
clc
rng(0);                                                                                 % We fix a seed to get my results replicable in case of randomness
warning('off');                                                                         % No warnings
options = optimset('Display','none');                                                   % Suppress function messages
set(0,'DefaultFigureWindowStyle','docked');                                             % Docking figures (more manageable)
format long                                                                             % More representation precision
FntNm = 'Times';
FntSz = 20;
formatData = 'dd/mm/yy'; 
IBDaycount = 3;

%% Calibration of Forwards and Discount factors

% Data Files
SX5Einputfile = '20230712SX5E.xlsx';
SPXinputfile  = '20230712SPX.xlsx';
N225inputfile = '20230712N225.xlsx';

% import data
if ispc()   
    [SX5EoptionData] = readExcelData(STOXX50Einputfile, formatData);
    [SPXoptionData]  = readExcelData(SPXinputfile, formatData);
    [N225optionData] = readExcelData(N225inputfile, formatData);
else    
    [SX5EoptionData] = readExcelDataMacOS(SX5Einputfile, "SX5EAll");
    [SPXoptionData]  = readExcelDataMacOS(SPXinputfile, "SPXAll");
    [N225optionData] = readExcelDataMacOS(N225inputfile, "N225All");
end

% Correlations
load('SPXSX5EHistCorr.mat')
% load('OverlapCorr.mat')


% names of the assets
AssetNames = { '.S&P500', '.STOXX50E', '.N225'};
SetDate    = SPXoptionData.settlement;

% Preprocessing
[SX5EDataTab, SPXDataTab, N225DataTab] = PreProc(SX5EoptionData, SPXoptionData, N225optionData);

% liquidity constraints
SPXDataTab  = liquidityConstraints(SPXDataTab);
SX5EDataTab = liquidityConstraints(SX5EDataTab);
N225DataTab = liquidityConstraints(N225DataTab);

% discount factors and forwards
[SPXDates, SPXDiscounts, SPXForwards]    = bootstrapMarketDiscountslm(SetDate, SPXDataTab, AssetNames{1,1});
[SX5EDates, SX5EDiscounts, SX5EForwards] = bootstrapMarketDiscountslm(SetDate, SX5EDataTab, AssetNames{1,2});
[N225Dates, N225Discounts, N225Forwards] = bootstrapMarketDiscountslm(SetDate, N225DataTab, AssetNames{1,3});

%% Calibration of on just one maturity

% find a maturity (reduce the dataset)
RequiredMat = 6/12; % maturity on which calibration is done
SetDate     = SPXoptionData.settlement;

% calibration of the model
q            = [0.5; 0.75; 1; 1.25; 1.5; 2];
SPXMarginalRMSE  = zeros(size(q,1),1);
SX5EMarginalRMSE = zeros(size(q,1),1);
SPXMarginalMAPE  = zeros(size(q,1),1);
SX5EMarginalMAPE = zeros(size(q,1),1);
CorrRMSE     = zeros(size(q,1),1);
x00 = [2, -4, 0.012];
TimeHorizons = [1; 5; 20; 60; 120];
CorrMatrix   = zeros(length(TimeHorizons), length(q)); 
% OverlapCorr  = [TimeHorizons, Correlations(:,2)];
Correlations   = [TimeHorizons, Correlations(:,2)];
BoundZero      = zeros(length(TimeHorizons), length(q));
BoundInf       = zeros(length(TimeHorizons), length(q));
SPXMargParams  = zeros(length(q), 3);
SX5EMargParams = zeros(length(q), 3);
rhoParam       = zeros(length(q), 1);
aParam         = zeros(length(q), 1);

for i=1:size(q,1)
    [MarginalParams, CommonParams, CalibError, CalibCorr, Bounds] = CalibrateSato1Mat(RequiredMat, SetDate, ...
        SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
        q(i), 0.5, x00, AssetNames, Correlations);
    
    rhoParam(i)  = CommonParams(1);
    aParam(i)    = CommonParams(2);
    
    SPXMargParams(i,:)  = MarginalParams(1,:);
    SX5EMargParams(i,:) = MarginalParams(2,:);
    
    CorrMatrix(:, i)    = CalibCorr;
    BoundZero(:, i)     = Bounds(1).*ones(size(TimeHorizons,1),1);
    BoundInf(:, i)      = Bounds(2).*ones(size(TimeHorizons,1),1);
    SPXMarginalRMSE(i)  = CalibError.Marginals.RMSE(1);
    SX5EMarginalRMSE(i) = CalibError.Marginals.RMSE(2);
    SPXMarginalMAPE(i)  = CalibError.Marginals.MAPE(1);
    SX5EMarginalMAPE(i) = CalibError.Marginals.MAPE(2);
    CorrRMSE(i)         = CalibError.Common.RMSE;
end

% Zmean  = sqrt(pi)*CommonParams(2);
% ZVar   = 0.5*sqrt(pi)*CommonParams(2);
% S1mean = sqrt(pi*MarginalParams(1,1));
% S1Var  = 0.5*sqr t(pi*MarginalParams(1,1)^3);
% S2mean = sqrt(pi*MarginalParams(2,1));
% S2Var  = 0.5*sqrt(pi*MarginalParams(2,1)^3);
% 
% rhoMIN = CommonParams(1)*sqrt(MarginalParams(1,1)*MarginalParams(2,1))*Zmean/sqrt(S1mean*S2mean);
% rhoMAX = (MarginalParams(1,1)*MarginalParams(2,1)*ZVar)/(sqrt(S1Var*S2Var));
% 
% rhoMINvec = rhoMIN.*ones(size(CorrMatrix));
% rhoMAXvec = rhoMAX.*ones(size(CorrMatrix));

CO1 = zeros(size(CorrMatrix,1), size(CorrMatrix,2)); % red
CO2 = ones(size(CorrMatrix,1), size(CorrMatrix,2)).*linspace(0.5,0.6,size(CorrMatrix,2)); % green
CO3 = ones(size(CorrMatrix,1), size(CorrMatrix,2)).*linspace(0,1,size(CorrMatrix,2)); % blue

figure
s1 = surf(q, TimeHorizons./365,CorrMatrix, 'DisplayName', 'CorrSurface');
s1.FaceColor = 'c';
s1.EdgeColor = 'r';
hold on
s2 = surf(q, TimeHorizons./365,BoundZero, CO2, 'DisplayName', 'BoundToZero');
s2.FaceColor = 'm';
s2.EdgeColor = 'b';
hold on
s3 = surf(q, TimeHorizons./365,BoundInf, CO3, 'DisplayName', 'BoundToInf');
s3.FaceColor = 'y';
s3.EdgeColor = 'g';
xlabel('q-Value', 'FontSize', FntSz, 'FontName',FntNm);
ylabel('TimeHorizon', 'FontSize', FntSz, 'FontName',FntNm)
hold off
legend

TimeHorizons = 1/365.*[1; 5; 20; 60; 120];

figure
ax = gca;
plot(TimeHorizons, 100*CorrMatrix(:,2), 'LineWidth', 2, 'DisplayName', 'Calibrated Correlation')
hold on
plot(TimeHorizons, 100*Correlations(:,2), 'LineWidth', 2, 'DisplayName', 'Historical Correlation')
grid on
legend
xlabel('Time Horizon')
ylabel('Correlation Values')
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';
ax.YLim = [50 100];


