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
% load('SPXSX5EHistCorr.mat')
% load('SPXN225HistCorr.mat')
load('SPXSX5EIncrCorr.mat');


% names of the assets
AssetNames1 = { '.S&P500', '.STOXX50E', '.N225'};
AssetNames = { '.S&P500', '.STOXX50E'};
SetDate    = SPXoptionData.settlement;

% Preprocessing
[SX5EDataTab, SPXDataTab, N225DataTab] = PreProc(SX5EoptionData, SPXoptionData, N225optionData);

% liquidity constraints
SPXDataTab  = liquidityConstraints(SPXDataTab);
SX5EDataTab = liquidityConstraints(SX5EDataTab);
N225DataTab = liquidityConstraints(N225DataTab);

% discount factors and forwards
[SPXDates, SPXDiscounts, SPXForwards]    = bootstrapMarketDiscountslm(SetDate, SPXDataTab, AssetNames1{1,1});
[SX5EDates, SX5EDiscounts, SX5EForwards] = bootstrapMarketDiscountslm(SetDate, SX5EDataTab, AssetNames1{1,2});
[N225Dates, N225Discounts, N225Forwards] = bootstrapMarketDiscountslm(SetDate, N225DataTab, AssetNames1{1,3});

%% Calibration on just one maturity

% find a maturity (reduce the dataset)
RequiredMat = 1.2/12; % maturity on which calibration is done
SetDate     = SPXoptionData.settlement;

% calibration of the model
% q = [0.5; 0.75; 1; 1.25];
q = 0.9;
x00 = [2, -4, 0.012];
MarginalParams = zeros(2, 3, length(q));
CommonParams   = zeros(length(q), 2);
qCal           = zeros(length(q), 1);
CorrError      = zeros(length(q), 1);
CorrMatrix     = zeros(5, length(q));
for i=1:length(q)
    [MarginalParams(:, :, i), CommonParams(i, :), Errors, CalibCorr, CorrBounds, qCal(i)] = CalibrateSato1Mat(RequiredMat, SetDate, ...
            SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
            q(i), 0.5, x00, AssetNames, Correlations);
%     [MarginalParams(:, :, i), CommonParams(i, :), Errors, CalibCorr, CorrBounds, qCal(i)] = CalibrateSato1Mat(RequiredMat, SetDate, ...
%             SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
%             q(i), 0.5, x00, AssetNames, Correlations);
    CorrError(i) = Errors.Common.RMSE;
    CorrMatrix(:, i) = CalibCorr;
end

BestCalCorr = load('SPXSX5EIncrCorr.mat');
BestCalCorr = BestCalCorr.Correlations(:, 2);
figure
ax = gca;
plot(Correlations(:,1), 100*CalibCorr, '*-', 'LineWidth', 2, 'DisplayName', 'CalibratedCorr')
hold on
plot(Correlations(:,1), 100*BestCalCorr, '*-', 'LineWidth', 2, 'DisplayName', 'HistCorr')
ylim([-100,100])
grid on
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.YAxis.TickLabelFormat = '%g%%';
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
xlabel('Time Horizon (Days)')
ylabel('Correlation Values')
legend
