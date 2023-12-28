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
formatData1 = 'dd/mm/yy'; 
formatData2 = 'mm/yy';
IBDaycount = 3;

%% Calibration of Forwards and Discount factors

% Data Files
% SX5EOptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230712SX5E.xlsx';
% SPXOptInputfile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230712SPX.xlsx';
% N225OptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230712N225.xlsx';
SX5EOptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230727SX5E.xlsx';
SPXOptInputfile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230727SPX.xlsx';
N225OptInputfile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/20230727N225.xlsx';

SPXSX5ECorrFile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/SPXSX5EPricesSheet.xlsm';
% SX5EN225CorrFile = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/SPXN225Prices.xlsx';
% SPXN225CorrFile  = '/Users/gabrielemasatti/Desktop/Magistrale/Master Thesis/Dataset/SPXN225PricesSheet.xlsm';

% Prices = xlsread(SPXSX5ECorrFile, "PricesFilteredTrue");
% Prices = xlsread(SX5EN225CorrFile, "PricesFilteredTrue");
% Prices = xlsread(SPXN225CorrFile, "PricesFilteredTrue");

% import data
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
SetDate    = SPXoptionData.settlement;

% Preprocessing
[SX5EDataTab, SPXDataTab, N225DataTab] = PreProc(SX5EoptionData, SPXoptionData, N225optionData);

% liquidity constraints and chosen dates
SPXDataTab  = liquidityConstraints(SPXDataTab);
SPXInitDate = '1-aug-2023';                         % modify this to take only portions of vol surface
SPXInitDate = datenum(SPXInitDate);
SPXEndDate  = '1-oct-2024';
SPXEndDate  = datenum(SPXEndDate);
SPXIndDates = find(SPXDataTab.MATURITIES >= SPXInitDate & SPXDataTab.MATURITIES <= SPXEndDate);
SPXDataTab  = SPXDataTab(SPXIndDates, :);
SX5EDataTab = liquidityConstraints(SX5EDataTab);
SX5EInitDate = '1-aug-2023';                         % modify this to take only portions of vol surface
SX5EInitDate = datenum(SX5EInitDate);
SX5EEndDate  = '1-oct-2024';
SX5EEndDate  = datenum(SX5EEndDate);
SX5EIndDates = find(SX5EDataTab.MATURITIES >= SX5EInitDate & SX5EDataTab.MATURITIES <= SX5EEndDate);
SX5EDataTab  = SX5EDataTab(SX5EIndDates, :);
N225DataTab = liquidityConstraints(N225DataTab);
N225InitDate = '1-aug-2023';                         % modify this to take only portions of vol surface
N225InitDate = datenum(N225InitDate);
N225EndDate  = '1-Jan-2024';
N225EndDate  = datenum(N225EndDate);
N225IndDates = find(N225DataTab.MATURITIES >= N225InitDate & N225DataTab.MATURITIES <= N225EndDate);
N225DataTab  = N225DataTab(N225IndDates, :);

% discount factors and forwards
[SPXDates, SPXDiscounts, SPXForwards]    = bootstrapMarketDiscountslm(SetDate, SPXDataTab, AssetNames{1,1});
[SX5EDates, SX5EDiscounts, SX5EForwards] = bootstrapMarketDiscountslm(SetDate, SX5EDataTab, AssetNames{1,2});
[N225Dates, N225Discounts, N225Forwards] = bootstrapMarketDiscountslm(SetDate, N225DataTab, AssetNames{1,3});
SPXRates  = - log(SPXDiscounts)./yearfrac(SetDate, SPXDates, IBDaycount);
SX5ERates = - log(SX5EDiscounts)./yearfrac(SetDate, SX5EDates, IBDaycount);
N225Rates = - log(N225Discounts)./yearfrac(SetDate, N225Dates, IBDaycount);

%% Calibration of NIG marginals on the vol surface
SetDate = SPXoptionData.settlement;
SPXS0   = 4607.41;
SX5ES0  = 4474.44;
%N225S0  = 4360.46;

% dividend yields
SPXdiv = SPXRates - ...
    (1./yearfrac(SetDate, SPXDates, IBDaycount).*log(SPXForwards/SPXS0));
SX5Ediv = SX5ERates - ...
    (1./yearfrac(SetDate, SX5EDates, IBDaycount).*log(SX5EForwards/SX5ES0));
% N225div = N225Rates - ...
%     (1./yearfrac(SetDate, N225Dates, IBDaycount).*log(N225Forwards/N225S0));


%% Historical Correlations
% Uncomment names of pairs used
AssetNames = {'.SPX'; '.SX5E'};
% AssetNames = {'.SPX'; '.N225'};

% TimeHorizons = [1; 5; 20; 60; 120; 240];
TimeHorizons = [ 5; 20; 60; 120; 240];

[Correlations, Returns] = AssetCorrelations(SPXSX5ECorrFile);
% [Correlations, Returns] = AssetCorrelations(SPXN225CorrFile);
Correlations = Correlations(2:end);             % without daily for N225
AnnualCorr = Correlations(end);
% Correlations = Correlations(1:end-1); 

figure();
ax = gca;
plot(TimeHorizons, 100*Correlations, '*-', 'DisplayName', 'Hist Correlation')
legend('FontName', FntNm, 'FontSize', FntSz)
grid on
title(strcat(AssetNames{1,1}, {' '}, AssetNames{2,1}, {' '}, 'Historical correlation vs TimeHorizon (in days)'), 'FontName', FntNm, 'FontSize', FntSz)
xlabel('TimeHorizon', 'FontName', FntNm, 'FontSize', FntSz)
ylabel('Correlation', 'FontName', FntNm, 'FontSize', FntSz)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.YAxis.TickLabelFormat = '%g%%';
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;

%% Calibration of Copula model on the vol surface
% Comment/Uncomment in the for loop to give in input the correct 
% datasets (referring to the required pair)

%%% OPTIMAL PARAMS %%%%
% q = 0.7;                    % for the pair .SPX .SX5E
q = 0.67;
% q1 = 0.65;                  % for the pair .SPX .N225 
% q2 = 0.4
% q = [0.4:0.02:1.4]';       % range checked (required ages to run)

% starting points
x00 = [0.1, -4, 1];       %[delta, beta, k]
% x00 = [0.1, -4, 1, 0.4];       %[delta, beta, k]
x00NIG = [2, 1, 0.012];       %[a, beta, delta]
CorrMatrix   = zeros(length(TimeHorizons), length(q)); 

% preallocations
MarginalParams1  = zeros(length(q), 3);
MarginalParams2  = zeros(length(q), 3);

% uncomment the proper ones to match the required pair
A1PricesMAPE    = zeros(length(q), length(SPXDates));
% A1PricesMAPE   = zeros(length(q), length(SX5EDates));
% A1PricesMAPE   = zeros(length(q), length(N225Dates));
A1VolsMAPE      = zeros(length(q), length(SPXDates));
% A1VolsMAPE     = zeros(length(q), length(SX5EDates));
% A1VolsMAPE     = zeros(length(q), length(N225Dates));

% A2PricesMAPE     = zeros(length(q), length(SPXDates));
A2PricesMAPE     = zeros(length(q), length(SX5EDates));
% A2PricesMAPE     = zeros(length(q), length(N225Dates));
% A2VolsMAPE      = zeros(length(q), length(SPXDates));
A2VolsMAPE     = zeros(length(q), length(SX5EDates));
% A2VolsMAPE     = zeros(length(q), length(N225Dates));

for i=1:length(q)
    [MarginalParams, CalibError, qcal] = CalibrateSatoNIGSurface(SetDate, ...
            SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
            q, 0.5, x00, AssetNames);
%     [MarginalParams, CalibError, qcal] = CalibrateSatoNIGSurface(SetDate, ...
%             SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
%              q, 0.5, x00, AssetNames);
    % [MarginalParams, CommonParams, CalibError, CalibCorr] = CalibrateSatoNIGSurface(SetDate, ...
    %         SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
    %         q, 0.5, x00, AssetNames, Correlations);

        %%% NIG Calibration %%%
%     [MarginalParamsNIG, CalibErrorNIG, ~, ~] = calibrateNIGcallput(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%             0.5, x00NIG, AssetNames);
%     [MarginalParamsNIG, CalibErrorNIG, ~, ~] = calibrateNIGcallput(SetDate, ...
%             SPXDiscounts, N225Discounts, SPXForwards, N225Forwards, SPXDataTab, N225DataTab, ...
%             0.5, x00NIG, AssetNames);
%     [MarginalParamsNIG, CalibErrorNIG, ~, ~] = calibrateNIGcallput(SetDate, ...
%             SPXDiscounts, SX5EDiscounts, SPXForwards, SX5EForwards, SPXDataTab, SX5EDataTab, ...
%             0.5, x00NIG, AssetNames);
    MarginalParams1(i, :)  = MarginalParams(1, :);
    MarginalParams2(i, :)  = MarginalParams(2, :);
    A1VolsMAPE(i, :)     = CalibError.volMAPE.Asset1';
    A2VolsMAPE(i, :)     = CalibError.volMAPE.Asset2';
    A1PricesMAPE(i, :)   = CalibError.PricesMAPE.Asset1';
    A2PricesMAPE(i, :)   = CalibError.PricesMAPE.Asset2';

end

% Graphical comparison with NIG

% .SPX
% figure
% ax = gca;
% plot(SPXDates, CalibError.PricesMAPE.Asset1, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIGSSD')
% hold on
% plot(SPXDates, CalibErrorNIG.PricesMAPE.Asset1, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.SPX MAPE Prices NIGSSD vs NIG', FontName=FntNm, FontSize=FntSz)
% 
% figure
% ax = gca;
% plot(SPXDates, CalibError.volMAPE.Asset1, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIGSSD')
% hold on
% plot(SPXDates, CalibErrorNIG.volMAPE.Asset1, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.SPX MAPE Volatilities NIGSSD vs NIG', FontName=FntNm, FontSize=FntSz)

% .SX5E
% figure
% ax = gca;
% plot(SX5EDates, CalibError.PricesMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIGSSD')
% hold on
% plot(SX5EDates, CalibErrorNIG.PricesMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.SX5E MAPE Prices NIGSSD vs NIG', FontName=FntNm, FontSize=FntSz)
% 
% figure
% ax = gca;
% plot(SX5EDates, CalibError.volMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIGSSD')
% hold on
% plot(SX5EDates, CalibErrorNIG.volMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.SX5E MAPE Volatilities NIGSSD vs NIG', FontName=FntNm, FontSize=FntSz)

% N225
% figure
% ax = gca;
% plot(N225Dates, CalibError.PricesMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIGSSD')
% hold on
% plot(N225Dates, CalibErrorNIG.PricesMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.N225 MAPE Prices NIGSSD vs NIG', FontName=FntNm, FontSize=FntSz)
% 
% figure
% ax = gca;
% plot(N225Dates, CalibError.volMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIGSSD')
% hold on
% plot(N225Dates, CalibErrorNIG.volMAPE.Asset2, '-*', 'LineWidth', 2, 'DisplayName', 'MAPE NIG')
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% legend
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';
% title('.N225 MAPE Volatilities NIGSSD vs NIG', FontName=FntNm, FontSize=FntSz)
% 
% % SPX
% figure
% ax = gca;
% bar(SPXDates, A1PricesMAPE, 'DisplayName', 'Prices MAPE')
% SPXTitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{1,1});
% legend
% title(SPXTitle, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTick = SPXDates;
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
% bar(SPXDates, A1VolsMAPE, 'DisplayName', 'Volatilities MAPE')
% SPXTitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{1,1});
% legend
% title(SPXTitle, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTick = SPXDates;
% datetick('x', formatData2, 'keeplimits')
% ax.XTickLabelRotation = 45;
% grid on
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';

% SX5E
% figure
% ax = gca;
% bar(SX5EDates, A2PricesMAPE, 'DisplayName', 'Prices MAPE')
% SX5ETitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
% legend
% title(SX5ETitle, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTickLabel  = datestr(SX5EDates, formatData2);
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
% bar(SX5EDates, A2VolsMAPE, 'DisplayName', 'Volatilities MAPE')
% SPXTitle = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
% legend
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
% figure
% ax = gca;
% bar(N225Dates, A2PricesMAPE, 'DisplayName', 'Prices MAPE')
% N225Title = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
% legend
% title(N225Title, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTickLabel  = datestr(N225Dates, formatData2);
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
% bar(N225Dates, A2VolsMAPE, 'DisplayName', 'Volatilities MAPE')
% N225Title = strcat('Marginal Calibration errors per Maturity', {' '}, AssetNames{2,1});
% legend
% title(N225Title, 'FontSize', FntSz, 'FontName', FntNm)
% ax.XTickLabel  = datestr(N225Dates, formatData2);
% ax.XTickLabelRotation = 45;
% grid on
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% ax.YAxis.TickLabelFormat = '%g%%';

% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %  

%% Test on the simulation

Nsim = 10000000;
TTM  = yearfrac(SetDate, SPXDates(4), IBDaycount);  %

% simulation with Gaussian NIG Copula
[X1sato, X2sato] = simulateNIGCOPULA(Nsim, TTM, 0, MarginalParams, 0.7, qcal);

% check of pricing of a call option through this method
K = 4600;
A1Fwd = SPXForwards(4).*exp(X1sato);
A2Fwd = SX5EForwards(4).*exp(X2sato);

A1EUCallDiscPayout = SPXDiscounts(4).*max((A1Fwd - K), 0);
A1EUCallPriceMC = mean(A1EUCallDiscPayout);

A2EUCallDiscPayout = SX5EDiscounts(4).*max((A2Fwd - K), 0);
A2EUCallPriceMC = mean(A2EUCallDiscPayout);

Params = FFTparameters(15, 0.0025, 1);
A1EUCallPriceFFT = CallPricesSatoNIGFFT(SPXForwards(4), SPXDiscounts(4), log(SPXForwards(4)/K), TTM, MarginalParams1, Params, 0.5, qcal);
A2EUCallPriceFFT = CallPricesSatoNIGFFT(SX5EForwards(4), SX5EDiscounts(4), log(SX5EForwards(4)/K), TTM, MarginalParams2, Params, 0.5, qcal);

% Simulation with increments
Ncoupons   = 4;                       % tot n of coupons
MonDates   = [SetDate; zeros(Ncoupons-1, 1)];    % last is the maturity (first is setdate) 
CouponFreq = 1;  % in months

X1path     = zeros(Nsim, Ncoupons);
X2path     = zeros(Nsim, Ncoupons);
F1path     = [SPXForwards(4)*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
meanF1path = [SPXForwards(4); zeros(Ncoupons-1,1)];
F2path     = [SX5EForwards(4)*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
meanF2path = [SX5EForwards(4); zeros(Ncoupons-1,1)];

for i = 2:Ncoupons
    MonDates(i) = busdate(datenum(addtodate(SetDate, CouponFreq*(i-1), 'month')));
    t = yearfrac(SetDate, MonDates(i), IBDaycount);
    s = yearfrac(SetDate, MonDates(i-1), IBDaycount);
    [X1path(:,i-1), X2path(:,i-1)] = simulateNIGCOPULA(Nsim, t, s, MarginalParams, 0.7, qcal);
    F1path(:, i)  = F1path(:, i-1).*exp(X1path(:,i-1));
    meanF1path(i) = mean(F1path(:, i));
    F2path(:, i)  = F2path(:, i-1).*exp(X2path(:,i-1));
    meanF2path(i) = mean(F2path(:, i));
end

[X1LastStep, X2LastStep] = simulateNIGCOPULA(Nsim, TTM, t, MarginalParams, 0.7, qcal);
A1FwdSteps = F1path(:, end).*exp(X1LastStep);
A2FwdSteps = F2path(:, end).*exp(X2LastStep);

A1EUCallDiscPayoutSteps = SPXDiscounts(4).*max((A1FwdSteps - K), 0);
A1EUCallPriceMCSteps = mean(A1EUCallDiscPayoutSteps);

A2EUCallDiscPayoutSteps = SX5EDiscounts(4).*max((A2FwdSteps - K), 0);
A2EUCallPriceMCSteps = mean(A2EUCallDiscPayoutSteps);

% MonDates(end) = SPXDates(6);
% t = yearfrac(SetDate, MonDates(end), IBDaycount);
% s = yearfrac(SetDate, MonDates(end-1), IBDaycount);
% [X1path(:,end-1), ~] = simulateNIGCOPULA(Nsim, t-s, MarginalParams, 0.1, qcal);
% F1path(:, end)  = F1path(:, end-1).*exp(X1path(:,end-1));

% Actually Pricing the product
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

rho = 0.72;
[Price] = MBRCWOFCOPULApricing(Nsim, ULsMktData, SetDate, Ncoupons, 3, 1, 0.02, 0.6, MarginalParams, rho, qcal);

