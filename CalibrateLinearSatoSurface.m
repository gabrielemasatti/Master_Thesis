function [MarginalParams, CommonParams, CalibError, CalibCorr] = CalibrateLinearSatoSurface(SetDate, ...
            A1Discounts, A2Discounts, A1Forwards, A2Forwards, A1DataTab, A2DataTab, ...
                q, alpha, Correlation, AssetNames)

%
% Function that calibrates the Linear Sato BB Model
% It gives both the marginal parameters and the common ones as output
%
% INPUT
% SetDate:      settlement date of the forward
% A1Discounts:  discount factors at relevant maturities for first asset
% A2Discounts:  discount factors at relevant maturities for second asset
% A1Forwards:   forward prices F(t0; T) at relevant maturities for first asset
% A2Forwards:   forward prices F(t0; T) at relevant maturities for second asset
% A1DataTab:    table with the data of options with first asset as UL
% A2DataTab:    table with the data of options with second asset as UL
% q:            self-similar parameter, at first fixed for the marginals,
% alpha:        alpha of the distribution considered (0.5 -> NIG)
% x00:          starting points
% AssetNames:   Names of the asset of which market data are calibrated
% Correlations: Historical correlations of the assets considered
%
% OUTPUT
% MarginalParams: [a_j, beta_j, delta_j] j=1,2 
% CommonParams:   q, a, rho common parameters driving the dependence
% CalibError:     struct with calibration errors
% CalibCorr:      calibrated values of correlation
%

% Preliminaries
FntSz      = 20;
FntNm      = 'Times';
IBDaycount = 3;

%% Preallocate variables
A1Maturities    = unique(A1DataTab.MATURITIES);
A2Maturities    = unique(A2DataTab.MATURITIES);
A1TTM           = yearfrac(SetDate, A1Maturities, IBDaycount);
A2TTM           = yearfrac(SetDate, A2Maturities, IBDaycount);
A1rates         = zeros(length(A1TTM), 1);
A2rates         = zeros(length(A2TTM), 1);
A1CallMoneyness = cell(length(A1TTM), 1);
A1PutMoneyness  = cell(length(A1TTM), 1);
A2CallMoneyness = cell(length(A2TTM), 1);
A2PutMoneyness  = cell(length(A2TTM), 1);
A1CallStrikes   = cell(length(A1TTM), 1);
A1PutStrikes    = cell(length(A1TTM), 1);
A2CallStrikes   = cell(length(A2TTM), 1);
A2PutStrikes    = cell(length(A2TTM), 1);
A1CallPrices    = cell(length(A1TTM), 1);
A1PutPrices     = cell(length(A1TTM), 1);
A2CallPrices    = cell(length(A2TTM), 1);
A2PutPrices     = cell(length(A2TTM), 1);

% Parameters for the FFT
Params = FFTparameters(15, 0.0025, 1);
A1nOptperMat = zeros(length(A1TTM), 1);
A2nOptperMat = zeros(length(A2TTM), 1);

% Model Prices Asset1
for i = 1:length(A1TTM)

    if (i==1)
        [A1OTMCall, A1OTMPut] = filterOTMoptions(A1TTM(i), SetDate, A1DataTab, A1Forwards, A1Discounts, A1Maturities);
        A1nOptperMat(i)       = length([A1OTMCall.STRIKES; A1OTMPut.STRIKES]);
        A1CallMoneyness{i}  = log(A1Forwards(i)./A1OTMCall.STRIKES);
        A1PutMoneyness{i}   = log(A1Forwards(i)./A1OTMPut.STRIKES);
        A1CallStrikes{i}    = A1OTMCall.STRIKES;
        A1PutStrikes{i}     = A1OTMPut.STRIKES;
        A1CallPrices{i}     = A1OTMCall.MID;
        A1PutPrices{i}      = A1OTMPut.MID;
        A1rates(i)          = -log(A1Discounts(i))/A1TTM(i);
        A1mktVolsCall       = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CallPrices{i});
        A1mktVolsPut        = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1PutPrices{i}, 'Class','put');
        A1MdlPricesCall     = @(cal) CallPricesLinearSatoFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal, Params, alpha, 1);
        A1MdlPricesPut      = @(cal) CallPricesLinearSatoFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal, Params, alpha, 1) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
        A1MdlImpvCall       = @(cal) blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesCall(cal)');
        A1MdlImpvPut        = @(cal) blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesPut(cal)', "Class", "put");
    
    else

        [A1OTMCallTemp, A1OTMPutTemp] = filterOTMoptions(A1TTM(i), SetDate, A1DataTab, A1Forwards, A1Discounts, A1Maturities);
        A1nOptperMat(i)     = length([A1OTMCallTemp.STRIKES; A1OTMPutTemp.STRIKES]);
        A1OTMCall           = [A1OTMCall; A1OTMCallTemp];
        A1OTMPut            = [A1OTMPut; A1OTMPutTemp];
        A1CallMoneyness{i}  = log(A1Forwards(i)./A1OTMCallTemp.STRIKES);
        A1PutMoneyness{i}   = log(A1Forwards(i)./A1OTMPutTemp.STRIKES);
        A1CallStrikes{i}    = A1OTMCallTemp.STRIKES;
        A1PutStrikes{i}     = A1OTMPutTemp.STRIKES;
        A1CallPrices{i}     = A1OTMCallTemp.MID;
        A1PutPrices{i}      = A1OTMPutTemp.MID;
        A1rates(i)          = -log(A1Discounts(i))/A1TTM(i);
        A1mktVolsCallTemp   = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CallPrices{i});
        A1mktVolsPutTemp    = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1PutPrices{i}, 'Class','put');
        A1mktVolsCall       = [A1mktVolsCall; A1mktVolsCallTemp];
        A1mktVolsPut        = [A1mktVolsPut; A1mktVolsPutTemp];
        A1MdlPricesCallTemp = @(cal) CallPricesLinearSatoFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal, Params, alpha, 1);
        A1MdlPricesPutTemp  = @(cal) CallPricesLinearSatoFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal, Params, alpha, 1) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
        A1MdlPricesCall     = @(cal) [A1MdlPricesCall(cal), A1MdlPricesCallTemp(cal)];
        A1MdlPricesPut      = @(cal) [A1MdlPricesPut(cal), A1MdlPricesPutTemp(cal)];
        A1MdlImpvCallTemp   = @(cal) blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesCallTemp(cal)');
        A1MdlImpvPutTemp    = @(cal) blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesPutTemp(cal)', "Class", "put");
        A1MdlImpvCall       = @(cal) [A1MdlImpvCall(cal); A1MdlImpvCallTemp(cal)];
        A1MdlImpvPut        = @(cal) [A1MdlImpvPut(cal); A1MdlImpvPutTemp(cal)];
    end
end

A1mktvols   = [A1mktVolsCall; A1mktVolsPut];
A1MdlPrices = @(cal) [A1MdlPricesCall(cal), A1MdlPricesPut(cal)];
A1MdlImpv   = @(cal) [A1MdlImpvCall(cal); A1MdlImpvPut(cal)];

% Model Prices Asset2
for i = 1:length(A2TTM)

    if (i==1)
        [A2OTMCall, A2OTMPut] = filterOTMoptions(A2TTM(i), SetDate, A2DataTab, A2Forwards, A2Discounts, A2Maturities);
        A2nOptperMat(i)       = length([A2OTMCall.STRIKES; A2OTMPut.STRIKES]);
        A2CallMoneyness{i}  = log(A2Forwards(i)./A2OTMCall.STRIKES);
        A2PutMoneyness{i}   = log(A2Forwards(i)./A2OTMPut.STRIKES);
        A2CallStrikes{i}    = A2OTMCall.STRIKES;
        A2PutStrikes{i}     = A2OTMPut.STRIKES;
        A2CallPrices{i}     = A2OTMCall.MID;
        A2PutPrices{i}      = A2OTMPut.MID;
        A2rates(i)          = -log(A2Discounts(i))/A2TTM(i);
        A2mktVolsCall       = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CallPrices{i});
        A2mktVolsPut        = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2PutPrices{i}, 'Class','put');
        A2MdlPricesCall     = @(cal) CallPricesLinearSatoFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal, Params, alpha, 2);
        A2MdlPricesPut      = @(cal) CallPricesLinearSatoFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal, Params, alpha, 2) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        A2MdlImpvCall       = @(cal) blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesCall(cal)');
        A2MdlImpvPut        = @(cal) blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesPut(cal)', "Class", "put");
    
    else
        [A2OTMCallTemp, A2OTMPutTemp] = filterOTMoptions(A2TTM(i), SetDate, A2DataTab, A2Forwards, A2Discounts, A2Maturities);
        A2nOptperMat(i)     = length([A2OTMCallTemp.STRIKES; A2OTMPutTemp.STRIKES]);
        A2OTMCall           = [A2OTMCall; A2OTMCallTemp];
        A2OTMPut            = [A2OTMPut; A2OTMPutTemp];
        A2CallMoneyness{i}  = log(A2Forwards(i)./A2OTMCallTemp.STRIKES);
        A2PutMoneyness{i}   = log(A2Forwards(i)./A2OTMPutTemp.STRIKES);
        A2CallStrikes{i}    = A2OTMCallTemp.STRIKES;
        A2PutStrikes{i}     = A2OTMPutTemp.STRIKES;
        A2CallPrices{i}     = A2OTMCallTemp.MID;
        A2PutPrices{i}      = A2OTMPutTemp.MID;
        A2rates(i)          = -log(A2Discounts(i))/A2TTM(i);
        A2mktVolsCallTemp   = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CallPrices{i});
        A2mktVolsPutTemp    = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2PutPrices{i}, 'Class','put');
        A2mktVolsCall       = [A2mktVolsCall; A2mktVolsCallTemp];
        A2mktVolsPut        = [A2mktVolsPut; A2mktVolsPutTemp];
        A2MdlPricesCallTemp = @(cal) CallPricesLinearSatoFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal, Params, alpha, 2);
        A2MdlPricesPutTemp  = @(cal) CallPricesLinearSatoFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal, Params, alpha, 2) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        A2MdlPricesCall     = @(cal) [A2MdlPricesCall(cal), A2MdlPricesCallTemp(cal)];
        A2MdlPricesPut      = @(cal) [A2MdlPricesPut(cal), A2MdlPricesPutTemp(cal)];
        A2MdlImpvCallTemp   = @(cal) blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesCallTemp(cal)');
        A2MdlImpvPutTemp    = @(cal) blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesPutTemp(cal)', "Class", "put");
        A2MdlImpvCall       = @(cal) [A2MdlImpvCall(cal); A2MdlImpvCallTemp(cal)];
        A2MdlImpvPut        = @(cal) [A2MdlImpvPut(cal); A2MdlImpvPutTemp(cal)];
    end
end

A2mktvols   = [A2mktVolsCall; A2mktVolsPut];
A2MdlPrices = @(cal) [A2MdlPricesCall(cal), A2MdlPricesPut(cal)];
A2MdlImpv   = @(cal) [A2MdlImpvCall(cal); A2MdlImpvPut(cal)];

% Mkt Prices
% Asset1
A1MktPricesCall = A1OTMCall.MID;
A1MktPricesPut  = A1OTMPut.MID;
A1MktPrices     = [A1MktPricesCall; A1MktPricesPut];
% Asset2
A2MktPricesCall = A2OTMCall.MID;
A2MktPricesPut  = A2OTMPut.MID;
A2MktPrices     = [A2MktPricesCall; A2MktPricesPut];

% N of Options
% Asset1
A1nCalls      = length(A1OTMCall.STRIKES);
A1nPuts       = length(A1OTMPut.STRIKES);
A1nOpt        = A1nCalls + A1nPuts;

% Asset2
A2nCalls      = length(A2OTMCall.STRIKES);
A2nPuts       = length(A2OTMPut.STRIKES);
A2nOpt        = A2nCalls + A2nPuts;

% calibration procedure
x00  = [0.1371, -0.1, 2.2548, 0.1278, -0.1, 1.1527, 0.6310]';
% [sigma, eta, k, q]
% LB1  = [0, -120, 0.00001]';
LB  = [0.1, -10, 0, 0.1, -10, 0, 0.4]';
UB  = [120, 100, 5, 120, 100, 5, 1.2]';
% LB2  = [0, -120, 0.00001]';
% LB2  = [0, -120, 0.00001, 0, -120, 0.00001]';
% UB2  = [120, 100, 5, 120, 100, 5]';
options = optimoptions('fmincon', 'Display','off');

function [c,ceq] = constraints (cal)  
  c   = [];
  ceq = [cal(1).^2/(cal(2).^2.*cal(3)) - cal(4).^2/(cal(5).^2.*cal(6))];
end

RMSE1 = @(cal) sqrt((1/A1nOpt) * sum((A1MdlPrices(cal)' - A1MktPrices).^2));

RMSE2 = @(cal) sqrt((1/A2nOpt) * sum((A2MdlPrices(cal)' - A2MktPrices).^2));

calibparams = fmincon(@(cal) 0.5.*(RMSE1(cal) + RMSE2(cal)), x00, [],[], [], [], LB, UB, @(cal) constraints(cal), options);

A1sigma = calibparams(1);
A1eta   = calibparams(2);
A1k     = calibparams(3);

A2sigma = calibparams(4);
A2eta   = calibparams(5);
A2k     = calibparams(6);

%% Visualize the results
A1AvgMAPEvol = zeros(length(A1TTM),1);
A2AvgMAPEvol = zeros(length(A2TTM),1);

A1AvgMAPEprices = zeros(length(A1TTM),1);
A2AvgMAPEprices = zeros(length(A2TTM),1);

for i = 1:length(A1TTM)
    
    A1CalPricesCall    = real(CallPricesLinearSatoFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), calibparams, Params, alpha, 1));
    A1CalPricesPut     = real(CallPricesLinearSatoFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), calibparams, Params, alpha, 1)) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})';
    A1mktvolsCall      = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CallPrices{i});
    A1CalVolCall       = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CalPricesCall');
    A1mktvolsPut       = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1PutPrices{i}, 'Class','put');
    A1CalVolPut        = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1CalPricesPut', 'Class', 'put');
    A1MktPricesTemp    = [A1CallPrices{i}; A1PutPrices{i}];
    A1CalPrices        = [A1CalPricesCall'; A1CalPricesPut'];
    A1mktvols          = [A1mktvolsCall; A1mktvolsPut];
    A1CalVol           = [A1CalVolCall; A1CalVolPut];
    Moneyness          = [-A1CallMoneyness{i}; -A1PutMoneyness{i}];
    A1AvgMAPEvol(i)    = 100.*(1/length(A1MktPricesTemp).*sum(abs(A1mktvols - A1CalVol)./A1mktvols));
    A1AvgMAPEprices(i) = 100.*(1/length(A1MktPricesTemp).*sum(abs(A1MktPricesTemp-A1CalPrices)./A1MktPricesTemp));

    figure
    ax = gca;
    plot(Moneyness, A1mktvols, '+', 'MarkerSize', 10 ,'DisplayName','MarketVols')
    hold on
    plot(Moneyness, A1CalVol, 'square' , 'MarkerSize', 10, 'DisplayName','CalibratedVols')
    ax.XAxis.FontSize = FntSz;
    ax.YAxis.FontSize = FntSz;
    ax.XAxis.FontName = FntNm;
    ax.YAxis.FontName = FntNm;
    legend
    ylim([0.05 0.3])
    titleName = strcat({'Linear Sato Smile for'}, {' '}, AssetNames{1,1}, {' '}, {'at'}, {' '}, {datestr(A1Maturities(i))});
    title(titleName, 'FontSize', FntSz, 'FontName', FntNm)
    xlabel('Moneyness', 'FontSize', FntSz, 'FontName', FntNm)
    ylabel('Implied Volatility', 'FontSize', FntSz, 'FontName', FntNm)
    grid on
    hold off

end

for i = 1:length(A2TTM)
    A2CalPricesCall    = real(CallPricesLinearSatoFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), calibparams, Params, alpha, 2));
    A2CalPricesPut     = real(CallPricesLinearSatoFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), calibparams, Params, alpha, 2)) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})';
    A2mktvolsCall      = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CallPrices{i});
    A2CalVolCall       = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CalPricesCall');
    A2mktvolsPut       = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2PutPrices{i}, 'Class','put');
    A2CalVolPut        = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2CalPricesPut', 'Class', 'put');
    A2MktPricesTemp    = [A2CallPrices{i}; A2PutPrices{i}];
    A2CalPrices        = [A2CalPricesCall'; A2CalPricesPut'];
    A2mktvols          = [A2mktvolsCall; A2mktvolsPut];
    A2CalVol           = [A2CalVolCall; A2CalVolPut];
    A2AvgMAPEvol(i)    = 100.*(1/length(A2MktPricesTemp).*sum(abs(A2mktvols - A2CalVol)./A2mktvols));
    A2AvgMAPEprices(i) = 100.*(1/length(A2MktPricesTemp).*sum(abs((A2MktPricesTemp-A2CalPrices))./A2MktPricesTemp));

    figure
    ax = gca;
    plot([-A2CallMoneyness{i}; -A2PutMoneyness{i}], A2mktvols, '+', 'MarkerSize', 10, 'DisplayName','MarketVols')
    hold on
    plot([-A2CallMoneyness{i}; -A2PutMoneyness{i}], A2CalVol, 'square' , 'MarkerSize', 10,'DisplayName','CalibratedVols')
    ylim([0.05 0.3])
    ax.XAxis.FontSize = FntSz;
    ax.YAxis.FontSize = FntSz;
    ax.XAxis.FontName = FntNm;
    ax.YAxis.FontName = FntNm;
    legend
    titleName = strcat({'Linear Sato Smile for '}, AssetNames{2,1}, {' '}, {'at'}, {' '}, {datestr(A2Maturities(i))});    
    title(titleName, 'FontSize', FntSz, 'FontName', FntNm)
    xlabel('Moneyness', 'FontSize', FntSz, 'FontName', FntNm)
    ylabel('Implied Volatility', 'FontSize', FntSz, 'FontName', FntNm)
    grid on
    hold off

end
A1CalPrices = A1MdlPrices(calibparams)';
A2CalPrices = A2MdlPrices(calibparams)';

%% Output on Marginals
% MarginalParams = struct;
% MarginalParams.Asset1 = [A1a; A1beta; A1delta];
% MarginalParams.Asset2 = [A2a; A2beta; A2delta];
MarginalParams = [A1sigma, A1eta, A1k; A2sigma, A2eta, A2k];
CalibError = struct;
CalibError.volRMSE.Asset1 = sqrt(1/A1nOpt.*sum((A1CalVol - A1mktvols).^2));
CalibError.volRMSE.Asset2 = sqrt(1/A2nOpt.*sum((A2CalVol - A2mktvols).^2));
CalibError.PricesRMSE.Asset1 = sqrt(1/A1nOpt.*sum((A1CalPrices - A1MktPrices).^2));
CalibError.PricesRMSE.Asset2 = sqrt(1/A2nOpt.*sum((A2CalPrices - A2MktPrices).^2));
CalibError.volMAPE.Asset1 = A1AvgMAPEvol;
CalibError.volMAPE.Asset2 = A2AvgMAPEvol;
CalibError.PricesMAPE.Asset1 = A1AvgMAPEprices;
CalibError.PricesMAPE.Asset2 = A2AvgMAPEprices;

%% Common parameters
c       = A1sigma^2/(A1eta^2*A1k);

% calibration of common parameter
CorrHat = @(k) (sqrt(A1k*A2k))/k;
d       = @(k) abs(CorrHat(k) - Correlation);
LB      = max(A1k, A2k);
UB      = 10;
k_z     = fmincon(@(k) d(k), 0.1, [], [], [], [], LB, UB, [], options);

eta_z   = - sqrt(1/(c*k_z));
sigma_z = 1;

A1k_1 = A1k*k_z/(k_z - A1k);
A2k_2 = A2k*k_z/(k_z - A2k);

a1 = A1eta/eta_z*(A1k_1/(k_z+A1k_1));
a2 = A2eta/eta_z*(A2k_2/(k_z+A2k_2));

A1sigma_1 = sqrt(k_z*a1^2/A1k_1);
A2sigma_2 = sqrt(k_z*a2^2/A2k_2);

A1eta_1 = k_z*a1*eta_z/A1k_1;
A2eta_2 = k_z*a2*eta_z/A2k_2;

CommonParams = [sigma_z, eta_z, k_z, calibparams(end)];

CalibCorr = CorrHat(k_z);

MarginalParams = [MarginalParams, [A1sigma_1; A2sigma_2], [A1eta_1; A2eta_2]...
    , [A1k_1; A2k_2], [a1; a2]];

end