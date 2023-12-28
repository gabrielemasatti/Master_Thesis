% Msc Thesis Masatti Gabriele
%
% Based on Patrizia Semeraro:
% Multivariate tempered stable additive subordination for financial models
%
% Supervisors: Michele Azzone, Roberto Baviera
% 
% -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- -*- 
%
% Academic Year 2022 - 2023
%
% Main for Calibration and Simulation of S-IG Subordinated Brownian Motion
% 

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
FntSz = 25;
LgndFntSize = 15;
formatData1 = 'dd/mm/yy'; 
formatData2 = 'mm/yy';
IBDaycount = 3;

%% Datasets
% Comment/Uncomment below proper input filtered files, paying attention 
% to dates and to path

% Option Files
% SX5EOptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230712SX5E.xlsx';
% SPXOptInputfile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230712SPX.xlsx';
% N225OptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230712N225.xlsx';
SX5EOptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230727SX5E.xlsx';
SPXOptInputfile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230727SPX.xlsx';
N225OptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230727N225.xlsx';

% Time series prices
% SPXSX5ECorrFile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/SPXSX5EPricesSheet.xlsm';
SPXN225CorrFile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/SPXN225PricesSheet.xlsm';

% import data (pay attention to computer settings, I only modified for mac)
if ispc()   
    [SX5EoptionData] = readExcelData(STOXX50Einputfile, formatData1);
    [SPXoptionData]  = readExcelData(SPXOptInputfile, formatData1);
    [N225optionData] = readExcelData(N225OptInputfile, formatData1);
else    
    [SX5EoptionData] = readExcelDataMacOS(SX5EOptInputfile, "SX5EAll");
    [SPXoptionData]  = readExcelDataMacOS(SPXOptInputfile, "SPXAll");
    [N225optionData] = readExcelDataMacOS(N225OptInputfile, "N225All");
end

% names of the assets
AssetNames = { '.S&P500', '.STOXX50E', '.N225'};
SetDate    = SPXoptionData.settlement;   % settlement, same for both assets

% Work with Tables with only necessary fields
[SX5EDataTab, SPXDataTab, N225DataTab] = PreProc(SX5EoptionData, SPXoptionData, N225optionData);

%% Liquidity Constraints and choosen Dates
% Modify InitDate and EndDate to keep just some parts of the surface
% Keep only portions that are consistent with both assets
% Initialised with just one dates

SPXDataTab  = liquidityConstraints(SPXDataTab);
SPXInitDate = '1-jul-2023';                         % modify this to take only portions of vol surface
SPXInitDate = datenum(SPXInitDate);
SPXEndDate  = '1-jan-2024';                         % modify this to take only portions of vol surface
SPXEndDate  = datenum(SPXEndDate);
SPXIndDates = find(SPXDataTab.MATURITIES >= SPXInitDate & SPXDataTab.MATURITIES <= SPXEndDate);
SPXDataTab  = SPXDataTab(SPXIndDates, :);

SX5EDataTab = liquidityConstraints(SX5EDataTab);
SX5EInitDate = '1-jul-2023';                        % modify this to take only portions of vol surface
SX5EInitDate = datenum(SX5EInitDate);
SX5EEndDate  = '1-oct-2024';                        % modify this to take only portions of vol surface
SX5EEndDate  = datenum(SX5EEndDate);
SX5EIndDates = find(SX5EDataTab.MATURITIES >= SX5EInitDate & SX5EDataTab.MATURITIES <= SX5EEndDate);
SX5EDataTab  = SX5EDataTab(SX5EIndDates, :);

N225DataTab  = liquidityConstraints(N225DataTab);
N225InitDate = '1-jun-2023';                        % modify this to take only portions of vol surface
N225InitDate = datenum(N225InitDate);   
N225EndDate  = '1-jan-2024';                        % modify this to take only portions of vol surface 
N225EndDate  = datenum(N225EndDate);
N225IndDates = find(N225DataTab.MATURITIES >= N225InitDate & N225DataTab.MATURITIES <= N225EndDate);
N225DataTab  = N225DataTab(N225IndDates, :);

% Find Discount Factors and Forwards @ relevant dates
[SPXDates, SPXDiscounts, SPXForwards]    = bootstrapMarketDiscountslm(SetDate, SPXDataTab, AssetNames{1,1});
[SX5EDates, SX5EDiscounts, SX5EForwards] = bootstrapMarketDiscountslm(SetDate, SX5EDataTab, AssetNames{1,2});
[N225Dates, N225Discounts, N225Forwards] = bootstrapMarketDiscountslm(SetDate, N225DataTab, AssetNames{1,3});

% Zero Rates
SPXRates  = - log(SPXDiscounts)./yearfrac(SetDate, SPXDates, IBDaycount);
SX5ERates = - log(SX5EDiscounts)./yearfrac(SetDate, SX5EDates, IBDaycount);
N225Rates = - log(N225Discounts)./yearfrac(SetDate, N225Dates, IBDaycount);

% Starting Values

% @ 12-07-2023
% SPXS0   = 4472.16;
% SX5ES0  = 4360.46;
% N225S0  = 31943.93;

% @ 27-07-2023
SPXS0   = 4607.41;
SX5ES0  = 4474.44;
N225S0  = 32891.16;

% dividend yields
SPXdiv = SPXRates - ...
    (1./yearfrac(SetDate, SPXDates, IBDaycount).*log(SPXForwards/SPXS0));
SX5Ediv = SX5ERates - ...
    (1./yearfrac(SetDate, SX5EDates, IBDaycount).*log(SX5EForwards/SX5ES0));
N225div = N225Rates - ...
    (1./yearfrac(SetDate, N225Dates, IBDaycount).*log(N225Forwards/N225S0));

% just quick check
% SPXS0 = SPXForwards./exp((SPXRates - SPXdiv).*yearfrac(SetDate, SPXDates, 3));

%% Historical Correlations
% Uncomment names of pairs used
% AssetNames = {'.SPX'; '.SX5E'};
AssetNames = {'.SPX'; '.N225'};

% TimeHorizons = [1; 5; 20; 60; 120; 240];
TimeHorizons = [5; 20; 60; 120; 240];                % without daily for N225

% [Correlations, Returns] = AssetCorrelations(SPXSX5ECorrFile);
[Correlations, Returns] = AssetCorrelations(SPXN225CorrFile);
Correlations = Correlations(2:end);             % without daily for N225

AnnualCorr = Correlations(end);
% Correlations = Correlations(1:end-1);

figure();
ax = gca;
plot(TimeHorizons, 100*Correlations, '*-', 'DisplayName', 'Hist Correlation')
legend('FontName', FntNm, 'FontSize', FntSz)
grid on
title(strcat(AssetNames{1,1}, {'/'}, AssetNames{2,1}, {' '}, 'Historical correlation vs TimeHorizon (in days)'), 'FontName', FntNm, 'FontSize', FntSz)
xlabel('TimeHorizon', 'FontName', FntNm, 'FontSize', FntSz)
ylabel('Correlation', 'FontName', FntNm, 'FontSize', FntSz)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.YAxis.TickLabelFormat = '%g%%';
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;

%% Calibration of S-IG model on the vol surface
% Comment/Uncomment in the for loop to give in input the correct 
% datasets (referring to the required pair)

%%% OPTIMAL PARAMS %%%%
% q = 0.94;                   % for the pair .SPX .SX5E
q = 1.0;                  % for the pair .SPX .N225 
% q = [0.4:0.02:1.4]';       % range checked (required ages to run)

% copula
qcopula =  0.7;             % for the pair .SPX .SX5E
% qcopula = 0.67;             % for the pair .SPX .N225 

% S-IG
x00 = [2, -4, 0.012];         %[a, beta, delta]
x00NIG = [2, 1, 0.012];       %[a, beta, delta]
CorrMatrix   = zeros(length(TimeHorizons), length(q)); 

% preallocations
MarginalParams1  = zeros(length(q), 3);
MarginalParams2  = zeros(length(q), 3);
CommonParameters = zeros(length(q), 3);
CorrRMSE         = zeros(length(q), 1);

% uncomment the proper ones to match the required pair
A1PricesMAPE    = zeros(length(q), length(SPXDates));
% A1PricesMAPE   = zeros(length(q), length(SX5EDates));
% A1PricesMAPE   = zeros(length(q), length(N225Dates));
A1VolsMAPE      = zeros(length(q), length(SPXDates));
% A1VolsMAPE     = zeros(length(q), length(SX5EDates));
% A1VolsMAPE     = zeros(length(q), length(N225Dates));

% A2PricesMAPE     = zeros(length(q), length(SPXDates));
% A2PricesMAPE     = zeros(length(q), length(SX5EDates));
A2PricesMAPE     = zeros(length(q), length(N225Dates));
% A2VolsMAPE      = zeros(length(q), length(SPXDates));
% A2VolsMAPE     = zeros(length(q), length(SX5EDates));
A2VolsMAPE     = zeros(length(q), length(N225Dates));

% designed for several q values (use only one otherwise will require AGES)
for i=1:length(q)
    
    %%% S-IG Calibration %%%
%     [MarginalParams, CommonParams, CalibError, CalibCorr, CorrBounds, SPXNoptPerMat, SX5ENoptPerMat] = CalibrateSatoSurface(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%             q, 0.5, x00, AssetNames, Correlations, TimeHorizons);
    [MarginalParams, CommonParams, CalibError, CalibCorr, CorrBounds, SPXNoptPerMat, N225NoptPerMat] = CalibrateSatoSurface(SetDate, ...
            SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
            q(i), 0.5, x00, AssetNames, Correlations, TimeHorizons);
%      [MarginalParams, CommonParams, CalibError, CalibCorr] = CalibrateSatoSurface(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%             q, 0.5, x00, AssetNames, Correlations, TimeHorizons);
    % [SPXeta, SPXk, SPXsigma, SPXMSE, SPXMAPE] = calibrateNIGcallput(SetDate, SPXDataTab, SPXForwards, SPXDiscounts, SPXDates);
    
    %%% Copula Calibration %%%
%     [MarginalParamsCop, CalibErrorCop, qcalCop] = CalibrateSatoNIGSurface(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%              qcopula, 0.5, x00copula, AssetNames);
%     [MarginalParamsCop, CalibErrorCop, qcalCop] = CalibrateSatoNIGSurface(SetDate, ...
%             SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
%              qcopula, 0.5, x00copula, AssetNames);

    %%% BB Model 1 %%%
%     [MarginalParamsBB, CommonParamsBB, CalibErrorBB, CalibCorrBB] = CalibrateLinearSatoSurface(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%                 q, 0.5, AnnualCorr, AssetNames);
    [MarginalParamsBB, CommonParamsBB, CalibErrorBB, CalibCorrBB] = CalibrateLinearSatoSurface(SetDate, ...
            SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
                q, 0.5, AnnualCorr, AssetNames);

    %%% BB Model 2 %%%
%     [MarginalParamsBB2, CommonParamsBB2, CalibErrorBB2, CalibCorrBB2] = CalibrateLinearSatoSurface2(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%                 q, 0.5, AnnualCorr, AssetNames);
    [MarginalParamsBB2, CommonParamsBB2, CalibErrorBB2, CalibCorrBB2] = CalibrateLinearSatoSurface2(SetDate, ...
            SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
                q, 0.5, AnnualCorr, AssetNames);

    %%% NIG Calibration %%%
%     [MarginalParamsNIG, CalibErrorNIG, ~, ~] = calibrateNIGcallput(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%             0.5, x00NIG, AssetNames);
    [MarginalParamsNIG, CalibErrorNIG, ~, ~] = calibrateNIGcallput(SetDate, ...
            SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
            0.5, x00NIG, AssetNames);
%     [MarginalParamsNIG, CalibErrorNIG, ~, ~] = calibrateNIGcallput(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%             0.5, x00NIG, AssetNames);

    MarginalParams1(i, :)  = MarginalParams(1, :);
    MarginalParams2(i, :)  = MarginalParams(2, :);
    CommonParameters(i, :) = CommonParams';
    A1VolsMAPE(i, :)     = CalibError.volMAPE.Asset1';
    A2VolsMAPE(i, :)     = CalibError.volMAPE.Asset2';
    A1PricesMAPE(i, :)   = CalibError.PricesMAPE.Asset1';
    A2PricesMAPE(i, :)   = CalibError.PricesMAPE.Asset2';
    CorrRMSE(i)          = CalibError.CorrRMSE;

end

% Graphical comparison

% .SPX
figure
ax = gca;
plot(SPXDates, CalibError.PricesMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE S-IG')
hold on
plot(SPXDates, CalibErrorNIG.PricesMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE NIG')
hold on
plot(SPXDates, CalibErrorBB.PricesMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato')
hold on
plot(SPXDates, CalibErrorBB2.PricesMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato D-IG')
datetick('x', formatData2, 'keeplimits')
ax.XTickLabelRotation = 45;
grid on
legend(fontname=FntNm, FontSize=LgndFntSize)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';
title('.SPX MAPE Prices S-IG, NIG and Linear Sato models', FontName=FntNm, FontSize=FntSz)

figure
ax = gca;
plot(SPXDates, CalibError.volMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE S-IG')
hold on
plot(SPXDates, CalibErrorNIG.volMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE NIG')
hold on
plot(SPXDates, CalibErrorBB.volMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato')
hold on
plot(SPXDates, CalibErrorBB2.volMAPE.Asset1, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato D-IG')
datetick('x', formatData2, 'keeplimits')
ax.XTickLabelRotation = 45;
grid on
legend(fontname=FntNm, FontSize=LgndFntSize)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';
title('.SPX MAPE Volatilities S-IG, NIG and Linear Sato models', FontName=FntNm, FontSize=FntSz)

% % .SX5E
% figure
% ax = gca;
% plot(SX5EDates, CalibError.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE S-IG')
% hold on
% plot(SX5EDates, CalibErrorNIG.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE NIG')
% hold on
% plot(SX5EDates, CalibErrorBB.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato')
% hold on
% plot(SX5EDates, CalibErrorBB2.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato D-IG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend(fontname=FntNm, FontSize=LgndFntSize)
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.SX5E MAPE Prices S-IG, NIG and Linear Sato models', FontName=FntNm, FontSize=FntSz)
% 
% figure
% ax = gca;
% plot(SX5EDates, CalibError.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE S-IG')
% hold on
% plot(SX5EDates, CalibErrorNIG.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE NIG')
% hold on
% plot(SX5EDates, CalibErrorBB.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato')
% hold on
% plot(SX5EDates, CalibErrorBB2.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato D-IG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend(fontname=FntNm, FontSize=LgndFntSize)
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.SX5E MAPE Volatilities S-IG, NIG and Linear Sato models', FontName=FntNm, FontSize=FntSz)

% N225
figure
ax = gca;
plot(N225Dates, CalibError.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE S-IG')
hold on
plot(N225Dates, CalibErrorNIG.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE NIG')
hold on
plot(N225Dates, CalibErrorBB.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato')
hold on
plot(N225Dates, CalibErrorBB2.PricesMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato D-IG')
datetick('x', formatData2, 'keeplimits')
ax.XTickLabelRotation = 45;
grid on
legend(fontname=FntNm, FontSize=LgndFntSize)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';
title('.N225 MAPE Prices S-IG, NIG and Linear Sato models', FontName=FntNm, FontSize=FntSz)

figure
ax = gca;
plot(SPXDates, CalibError.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE S-IG')
hold on
plot(SPXDates, CalibErrorNIG.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE NIG')
hold on
plot(SPXDates, CalibErrorBB.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato')
hold on
plot(SPXDates, CalibErrorBB2.volMAPE.Asset2, '-*', 'LineWidth', 1.3, 'DisplayName', 'MAPE Linear Sato D-IG')
datetick('x', formatData2, 'keeplimits')
ax.XTickLabelRotation = 45;
grid on
legend(fontname=FntNm, FontSize=LgndFntSize)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';
title('.N225 MAPE Volatilities S-IG, NIG and Linear Sato models', FontName=FntNm, FontSize=FntSz)

%% Marginal Calibration errors per Mat
% Uncomment and comment plots to get the required pair results

% SPX
figure
ax = gca;
b = bar(SPXDates, CalibError.PricesMAPE.Asset1, 'DisplayName', 'Prices MAPE');
b.FaceColor = '#6D899E';
b.EdgeColor = '#6D899E';
SPXTitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{1,1});
legend(fontname=FntNm, FontSize=LgndFntSize)
title(SPXTitle, 'FontSize', FntSz, 'FontName', FntNm)
ax.XTick = SPXDates;
datetick('x', formatData2, 'keeplimits')
ax.XTickLabelRotation = 45;
grid on
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';

figure
ax = gca;
b = bar(SPXDates, CalibError.volMAPE.Asset1, 'DisplayName', 'Volatilities MAPE');
b.FaceColor = '#6D899E';
b.EdgeColor = '#6D899E';
SPXTitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{1,1});
legend(fontname=FntNm, FontSize=LgndFntSize)
title(SPXTitle, 'FontSize', FntSz, 'FontName', FntNm)
ax.XTick = SPXDates;
datetick('x', formatData2, 'keeplimits')
ax.XTickLabelRotation = 45;
grid on
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';

% % SX5E
% figure
% ax = gca;
% b = bar(SX5EDates, CalibError.PricesMAPE.Asset2, 'DisplayName', 'Prices MAPE');
% b.FaceColor = '#6D899E';
% b.EdgeColor = '#6D899E';
% SX5ETitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
% legend(fontname=FntNm, FontSize=LgndFntSize)
% title(SX5ETitle, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTick = SX5EDates;
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% 
% figure
% ax = gca;
% b = bar(SX5EDates, CalibError.volMAPE.Asset2, 'DisplayName', 'Volatilities MAPE');
% b.FaceColor = '#6D899E';
% b.EdgeColor = '#6D899E';
% SPXTitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
% legend(fontname=FntNm, FontSize=LgndFntSize)
% title(SX5ETitle, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTick = SX5EDates;
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';

% N225
figure
ax = gca;
b = bar(N225Dates, CalibError.PricesMAPE.Asset2, 'DisplayName', 'PricesMAPE');
b.FaceColor = '#6D899E';
b.EdgeColor = '#6D899E';
N225Title = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
legend(fontname=FntNm, FontSize=LgndFntSize)
title(N225Title, 'FontSize', FntSz, 'FontName', FntNm)
ax.XTickLabel  = datestr(N225Dates, formatData2);
ax.XTickLabelRotation = 45;
grid on
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';

figure
ax = gca;
b = bar(N225Dates, CalibError.volMAPE.Asset2, 'DisplayName', 'Volatilities MAPE');
b.FaceColor = '#6D899E';
b.EdgeColor = '#6D899E';
N225Title = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
legend(fontname=FntNm, FontSize=LgndFntSize)
title(N225Title, 'FontSize', FntSz, 'FontName', FntNm)
ax.XTickLabel  = datestr(N225Dates, formatData2);
ax.XTickLabelRotation = 45;
grid on
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';

%% Correlation fit
figure
ax = gca;
plot(TimeHorizons, 100*Correlations, '*-', 'DisplayName', 'Historical Correlation')
hold on
plot(TimeHorizons, 100*CalibCorr, '*-', 'DisplayName', 'Calibrated Correlation')
legend('FontName',FntNm, 'FontSize', FntSz)
grid on
title('Reproduction of the Historical Correlation', FontName=FntNm, FontSize=FntSz)
xlabel('TimeHorizons')
ylabel('Correlation Values')
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.YAxis.TickLabelFormat = '%g%%';
ax.YLim = [20 100];

% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %  

%% Check of the Pricing of a single Call
Nsim = 1e7; % n of simulations
TTM  = yearfrac(SetDate, SPXDates(11), IBDaycount);  % 6 months

% simulation with S-IG with Subordinator S (one-step)
[X1sato_Cont, X2sato_Cont] = simulateSATONIG(Nsim, TTM, MarginalParams, CommonParameters, 0);

% check pricing of a call option through this method
K = 4600;

% Forward Prices
A1Fwd_Cont = SPXForwards(11).*exp(X1sato_Cont);
A2Fwd_Cont = SX5EForwards(9).*exp(X2sato_Cont);

SPXEUCallDiscPayout_Cont = SPXDiscounts(11)*max((A1Fwd_Cont - K), 0);
SX5EEUCallDiscPayout_Cont = SX5EDiscounts(9)*max((A2Fwd_Cont - K), 0);

% Price of the two call options with MC
SPXEUCallPriceMC_Cont = mean(SPXEUCallDiscPayout_Cont);
SX5EUCallPriceMC_Cont = mean(SX5EEUCallDiscPayout_Cont);

% Price of the two call options with close formula
Params = FFTparameters(15, 0.0025, 1);
SPXEUCallPriceFFT  = CallPricesSatoFFT(SPXForwards(11), SPXDiscounts(11), log(SPXForwards(11)/K), TTM, MarginalParams1, Params, 0.5, CommonParameters(3));
SX5EEUCallPriceFFT = CallPricesSatoFFT(SX5EForwards(9), SX5EDiscounts(9), log(SX5EForwards(9)/K), TTM, MarginalParams2, Params, 0.5, CommonParameters(3));

% Simulate with FFT with no a, rho
X1 = simulateSATONIGFFT(Nsim, 12, 1, TTM, 0, MarginalParams1(1), MarginalParams1(2), MarginalParams1(3), CommonParameters(3));
X2 = simulateSATONIGFFT(Nsim, 12, 1, TTM, 0, MarginalParams2(1), MarginalParams2(2), MarginalParams2(3), CommonParameters(3));
A1Fwd = SPXForwards(11).*exp(X1);
A2Fwd = SX5EForwards(9).*exp(X2);
SPXEUCallDiscPayout  = SPXDiscounts(11)*max((A1Fwd - K), 0);
SX5EEUCallDiscPayout = SX5EDiscounts(9)*max((A2Fwd - K), 0);
SPXEUCallPriceMC  = mean(SPXEUCallDiscPayout);
SX5EEUCallPriceMC = mean(SX5EEUCallDiscPayout);

% Simulate with FFT with a and rho
%%% CORRECT PARAMS %%%
%%% M = 14 - b = 0.01 - omega = 0.9 - g2 = -3; %%%

[f1check, f2check] = simulateLogFwd(Nsim, 14, 1, TTM, 0, MarginalParams, CommonParameters);
A1Fwdcheck = SPXForwards(11).*exp(f1check);
A2Fwdcheck = SX5EForwards(9).*exp(f2check);
SPXEUCallDiscPayoutcheck  = SPXDiscounts(11)*max((A1Fwdcheck - K), 0);
SX5EEUCallDiscPayoutcheck = SX5EDiscounts(9)*max((A2Fwdcheck - K), 0);
SPXEUCallPriceMCcheck  = mean(SPXEUCallDiscPayoutcheck);
SX5EEUCallPriceMCcheck = mean(SX5EEUCallDiscPayoutcheck);

% Try to simulate a path step by step and check the path equality holds
Ncoupons   = 4;                       % tot n of coupons
MonDates   = [SetDate; zeros(Ncoupons-1, 1)];    % last is the maturity (first is setdate) 
CouponFreq = 3;  % in months

X1path     = zeros(Nsim, Ncoupons);
X2path     = zeros(Nsim, Ncoupons);
F1path     = [SPXForwards(11)*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
F2path     = [SX5EForwards(9)*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
meanF1disjoint     = [mean(SPXForwards(11)*ones(Nsim,1)); zeros(Ncoupons-1, 1)];

for i = 2:Ncoupons
    MonDates(i) = busdate(datenum(addtodate(SetDate, CouponFreq*(i-1), 'month')));
    t = yearfrac(SetDate, MonDates(i), IBDaycount);
    s = yearfrac(SetDate, MonDates(i-1), IBDaycount);
    X1path(:,i-1) = simulateSATONIGFFT(Nsim, 15, 1, t, s, MarginalParams1(1), MarginalParams1(2), MarginalParams1(3), CommonParameters(3));
    F1path(:, i)  = F1path(:, i-1).*exp(X1path(:,i-1));
    meanF1disjoint(i) = mean(F1path(:, i));
end

% simulate using IG subordinators (to exploit dependences)
%%% CORRECT PARAMS %%%
%%% M = 14 - b = 0.001 - omega = 0.9 - g2 = -1; %%%
Ncoupons   = 4;                       % tot n of coupons
MonDates   = [SetDate; zeros(Ncoupons-1, 1)];    % last is the maturity (first is setdate) 
CouponFreq = 3;  % in months

f1     = zeros(Nsim, Ncoupons);
f2     = zeros(Nsim, Ncoupons);
F1path = [SPXForwards(11)*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
F2path = [SX5EForwards(9)*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
meanF1 = [mean(SPXForwards(11)*ones(Nsim,1)); zeros(Ncoupons-1, 1)];
meanF2 = [mean(SX5EForwards(9)*ones(Nsim,1)); zeros(Ncoupons-1, 1)];

for i = 2:Ncoupons
    MonDates(i) = busdate(datenum(addtodate(SetDate, CouponFreq*(i-1), 'month')));
    t = yearfrac(SetDate, MonDates(i), IBDaycount);
    s = yearfrac(SetDate, MonDates(i-1), IBDaycount);
    [f1(:,i-1), f2(:,i-1)] = simulateLogFwd(Nsim, 14, 1, t, s, MarginalParams, CommonParameters);
    F1path(:, i)  = F1path(:, i-1).*exp(f1(:,i-1));
    F2path(:, i)  = F2path(:, i-1).*exp(f2(:,i-1));
    meanF1(i) = mean(F1path(:, i));
    meanF2(i) = mean(F2path(:, i));
end

% plot of a distribution of an increment
figure
ax = gca;
[pdfy, pdfx] = ksdensity(f1(:,end-1));
plot(pdfx, pdfy, 'LineWidth', 2, 'DisplayName', 'pdf')
legend
grid on 
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;
ax.XLim = [-0.4 0.4];
title('pdf of a forward increment', FontSize=FntSz, FontName=FntNm)

[f1LastStep, f2LastStep] = simulateLogFwd(Nsim, 15, 1, TTM, t, MarginalParams, CommonParameters);
A1FwdSteps = F1path(:, end).*exp(f1LastStep);
A2FwdSteps = F2path(:, end).*exp(f2LastStep);

A1EUCallDiscPayoutSteps = SPXDiscounts(11).*max((A1FwdSteps - K), 0);
A1EUCallPriceMCSteps = mean(A1EUCallDiscPayoutSteps);

A2EUCallDiscPayoutSteps = SX5EDiscounts(9).*max((A2FwdSteps - K), 0);
A2EUCallPriceMCSteps = mean(A2EUCallDiscPayoutSteps);

% check it coincides with the forward of the 1-step simulation
TTMcheck = yearfrac(SetDate, SPXDates(11), IBDaycount);  % 1 year
[X1sato_Contcheck, X2sato_Contcheck] = simulateSATONIG(Nsim, TTMcheck, MarginalParams, CommonParameters, 0);
F1check = SPXForwards(11).*exp(X1sato_Contcheck);

disp(mean(F1path(:, end)) - mean(F1check))

% Actually Pricing the Product
ULsMktData = struct;
ULsMktData.Asset1.S0 = SPXS0;
ULsMktData.Asset1.Dates = SPXDates;
ULsMktData.Asset1.Forwards  = SPXForwards;
ULsMktData.Asset1.Discounts = SPXDiscounts;
ULsMktData.Asset1.Rates     = SPXRates;
ULsMktData.Asset1.Div       = SPXdiv;
ULsMktData.Asset2.S0 = SX5ES0;
ULsMktData.Asset2.Dates = SX5EDates;
ULsMktData.Asset2.Forwards  = SX5EForwards;
ULsMktData.Asset2.Discounts = SX5EDiscounts;
ULsMktData.Asset2.Rates     = SX5ERates;
ULsMktData.Asset2.Div       = SX5Ediv;

[Price] = MBRCWOFpricing(Nsim, ULsMktData, SetDate, Ncoupons, 3, 1, 0.02, 0.6, MarginalParams, CommonParameters);

% rhoCopula = 0.855;
% % [PriceCopula, CorrSim] = MBRCWOFCOPULApricing(Nsim, ULsMktData, SetDate, Ncoupons, 3, 1, 0.02, 0.6, MarginalParamsCop, rhoCopula, qcalCop);
% 
% % check on the correlation structure
% [X1satoCop, X2satoCop] = simulateNIGCOPULA(Nsim, 1, 0, MarginalParamsCop, rhoCopula, qcalCop);
% A1FwdCop = SPXForwards(11).*exp(X1satoCop);
% A2FwdCop = SX5EForwards(9).*exp(X2satoCop);
% CorrCopulaEmp = corr(X1satoCop, X2satoCop);
% 
% % Trying with several rho values
% RhoValues = [0.2:0.1:0.8];
% aValues   = linspace(0.05, 0.1555, 7);
% PriceGridrho = zeros(length(RhoValues), 1);
% PriceGrida   = zeros(length(RhoValues), 1);
% 
% PriceGridCopula = zeros(length(RhoValues), 1);
% 
% for i=1:length(RhoValues)
%     [PriceGridrho(i)] = MBRCWOFpricing(Nsim, ULsMktData, SetDate, Ncoupons, 3, 1, 0.02, 0.6, MarginalParams, [RhoValues(i), CommonParameters(2), CommonParameters(3)]);
%     % [PriceGrida(i)]   = MBRCWOFpricing(Nsim, ULsMktData, SetDate, Ncoupons, 3, 1, 0.02, 0.6, MarginalParams, [CommonParameters(1), aValues(i), CommonParameters(3)]);
%     [PriceGridCopula(i)] = MBRCWOFCOPULApricing(Nsim, ULsMktData, SetDate, Ncoupons, 3, 1, 0.02, 0.6, MarginalParamsCop, RhoValues(i), qcalCop);
% end
% 
% figure
% ax = gca;
% plot(RhoValues, 100*PriceGridrho, '-*', 'LineWidth', 1.5, 'DisplayName', 'Prices S-IG Model')
% hold on
% plot(RhoValues, 100*PriceGridCopula, '-*', 'LineWidth', 1.5, 'DisplayName', 'Prices Copula Model')
% grid on
% legend(fontsize=FntSz, FontName=FntNm)
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% xlabel('rho value')
% ylabel('MBRCWOF Price')