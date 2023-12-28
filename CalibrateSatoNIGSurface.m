function [MarginalParams, CalibError, qcal] = CalibrateSatoNIGSurface(SetDate, ...
            A1Discounts, A2Discounts, A1Forwards, A2Forwards, A1DataTab, A2DataTab, ...
                  q, alpha, x00, AssetNames)
%
% Function that calibrates the NIG Sato model on the whole volatility surface.
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
% alpha:        alpha of the ETaS distribution considered (0.5 -> NIG)
% x00:          starting points
% AssetNames:   Names of the asset of which market data are callibrated
%
% OUTPUT
% MarginalParams: [delta_j, beta_j, k_j] j=1,2 as struct
% CalibError:     struct with calibration errors
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

%% Calibration of marginal params

% Parameters for the FFT
Params = FFTparameters(15, 0.0025, 1);

% Model Prices Asset1
for i = 1:length(A1TTM)

    if (i==1)
        [A1OTMCall, A1OTMPut] = filterOTMoptions(A1TTM(i), SetDate, A1DataTab, A1Forwards, A1Discounts, A1Maturities);
        A1CallMoneyness{i}  = log(A1Forwards(i)./A1OTMCall.STRIKES);
        A1PutMoneyness{i}   = log(A1Forwards(i)./A1OTMPut.STRIKES);
        A1CallStrikes{i}    = A1OTMCall.STRIKES;
        A1PutStrikes{i}     = A1OTMPut.STRIKES;
        A1CallPrices{i}     = A1OTMCall.MID;
        A1PutPrices{i}      = A1OTMPut.MID;
        A1rates(i)          = -log(A1Discounts(i))/A1TTM(i);
        A1mktVolsCall       = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CallPrices{i});
        A1mktVolsPut        = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1PutPrices{i}, 'Class','put');
        A1MdlPricesCall     = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal, Params, alpha, q);
        % A1MdlPricesCall     = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal, Params, alpha);
        A1MdlPricesPut      = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal, Params, alpha, q) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
        % A1MdlPricesPut      = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal, Params, alpha) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
        A1MdlImpvCall       = @(cal) blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesCall(cal)');
        A1MdlImpvPut        = @(cal) blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesPut(cal)', "Class", "put");
    
    else

        [A1OTMCallTemp, A1OTMPutTemp] = filterOTMoptions(A1TTM(i), SetDate, A1DataTab, A1Forwards, A1Discounts, A1Maturities);
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
        A1MdlPricesCallTemp = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal, Params, alpha, q);
        % A1MdlPricesCallTemp = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal, Params, alpha);
        A1MdlPricesPutTemp  = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal, Params, alpha, q) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
        % A1MdlPricesPutTemp  = @(cal) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal, Params, alpha) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
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
        A2CallMoneyness{i} = log(A2Forwards(i)./A2OTMCall.STRIKES);
        A2PutMoneyness{i}   = log(A2Forwards(i)./A2OTMPut.STRIKES);
        A2CallStrikes{i}    = A2OTMCall.STRIKES;
        A2PutStrikes{i}     = A2OTMPut.STRIKES;
        A2CallPrices{i}     = A2OTMCall.MID;
        A2PutPrices{i}      = A2OTMPut.MID;
        A2rates(i)          = -log(A2Discounts(i))/A2TTM(i);
        A2mktVolsCall       = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CallPrices{i});
        A2mktVolsPut        = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2PutPrices{i}, 'Class','put');
        A2MdlPricesCall     = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal, Params, alpha, q);
        % A2MdlPricesCall     = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal, Params, alpha);
        A2MdlPricesPut      = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal, Params, alpha, q) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        % A2MdlPricesPut      = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal, Params, alpha) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        A2MdlImpvCall       = @(cal) blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesCall(cal)');
        A2MdlImpvPut        = @(cal) blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesPut(cal)', "Class", "put");
    else
        [A2OTMCallTemp, A2OTMPutTemp] = filterOTMoptions(A2TTM(i), SetDate, A2DataTab, A2Forwards, A2Discounts, A2Maturities);
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
        A2MdlPricesCallTemp = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal, Params, alpha, q);
        % A2MdlPricesCallTemp = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal, Params, alpha);
        A2MdlPricesPutTemp  = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal, Params, alpha, q) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        % A2MdlPricesPutTemp  = @(cal) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal, Params, alpha) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
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
% LB1  = [0, -120, 0.00001]';
LB1  = [0, -10, 2]';
% LB1  = [0, -10, 2, 0]';
UB1  = [1, 5, 20]';
% UB1  = [1, 5, 20, 1.5]';
% LB2  = [0, -120, 0.00001]';
LB2  = [0, -120, 0.01]';
UB2  = [1, 5, 20]';
% LB2  = [0, -120, 0.01, 0]';
% UB2  = [1, 5, 20, 1.5]';
options = optimoptions('fmincon', 'Display','off');

aux1    = @(cal) sqrt((1/A1nOpt) * sum((A1MdlPrices(cal)' - A1MktPrices).^2));
cal1    = fmincon(@(cal) aux1(cal), x00, [],[], [], [], LB1, UB1, [], options);
% aux1    = @(cal) sqrt((1/A1nOpt) * sum((A1MdlImpv(cal) - A1mktvols).^2));
% cal1    = fmincon(@(cal) aux1(cal), x00, [],[], [], [], LB1, UB1, [], options);

aux2    = @(cal) sqrt((1/A2nOpt) * sum((A2MdlPrices(cal)' - A2MktPrices).^2));
cal2    = fmincon(@(cal) aux2(cal), x00, [],[], [], [], LB2, UB2, [], options);
% aux2    = @(cal) sqrt((1/A2nOpt) * sum((A2MdlImpv(cal) - A2mktvols).^2));
% cal2    = fmincon(@(cal) aux2(cal), x00, [],[], [], [], LB2, UB2, [], options);

% Asset 1 Params
A1delta = cal1(1);
A1theta = cal1(2);
A1k     = cal1(3);
% A1q     = cal1(4);

% Asset 2 Params
A2delta = cal2(1);
A2theta = cal2(2);
A2k     = cal2(3);
% A2q     = cal1(4);

%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-%-5-

%% Calibration of q

% Model Prices Asset1
for i = 1:length(A1TTM)

    if (i==1)

        A1MdlPricesCall = @(p) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal1, Params, alpha, p);
        A1MdlPricesPut  = @(p) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal1, Params, alpha, p) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})';
        A1MdlImpvCall   = @(p) blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesCall(p)');
        A1MdlImpvPut    = @(p) blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesPut(p)', "Class", "put");

    else

        A1MdlPricesCallTemp = @(p) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal1, Params, alpha, p);
        A1MdlPricesPutTemp  = @(p) CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal1, Params, alpha, p) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})'; 
        A1MdlPricesCall     = @(p) [A1MdlPricesCall(p), A1MdlPricesCallTemp(p)];
        A1MdlPricesPut      = @(p) [A1MdlPricesPut(p), A1MdlPricesPutTemp(p)];
        A1MdlImpvCallTemp   = @(p) blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesCallTemp(p)');
        A1MdlImpvPutTemp    = @(p) blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1MdlPricesPutTemp(p)', "Class", "put");
        A1MdlImpvCall       = @(p) [A1MdlImpvCall(p); A1MdlImpvCallTemp(p)];
        A1MdlImpvPut        = @(p) [A1MdlImpvPut(p); A1MdlImpvPutTemp(p)];
    
    end

end
A1MdlPrices = @(p) [A1MdlPricesCall(p), A1MdlPricesPut(p)];
A1RMSE      = @(p) sqrt(1/A1nOpt * sum((A1MdlPrices(p)' - A1MktPrices).^2));

A1MdlImpv   = @(p) [A1MdlImpvCall(p); A1MdlImpvPut(p)];
% A1RMSE      = @(p) sqrt(1/A1nOpt * sum((A1MdlImpv(p) - A1mktvols).^2));

% Model Prices Asset2
for i = 1:length(A2TTM)
    
    if (i==1)
        A2MdlPricesCall = @(p) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal2, Params, alpha, p);
        A2MdlPricesPut  = @(p) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal2, Params, alpha, p) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        A2MdlImpvCall   = @(p) blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesCall(p)');
        A2MdlImpvPut    = @(p) blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesPut(p)', "Class", "put");

    else

        A2MdlPricesCallTemp = @(p) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal2, Params, alpha, p);
        A2MdlPricesPutTemp  = @(p) CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal2, Params, alpha, p) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})'; 
        A2MdlPricesCall     = @(p) [A2MdlPricesCall(p), A2MdlPricesCallTemp(p)];
        A2MdlPricesPut      = @(p) [A2MdlPricesPut(p), A2MdlPricesPutTemp(p)];
        A2MdlImpvCallTemp   = @(p) blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesCallTemp(p)');
        A2MdlImpvPutTemp    = @(p) blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2MdlPricesPutTemp(p)', "Class", "put");
        A2MdlImpvCall       = @(p) [A2MdlImpvCall(p); A2MdlImpvCallTemp(p)];
        A2MdlImpvPut        = @(p) [A2MdlImpvPut(p); A2MdlImpvPutTemp(p)];
    
    end

end
A2MdlPrices     = @(p) [A2MdlPricesCall(p), A2MdlPricesPut(p)];
A2RMSE          = @(p) sqrt(1/A2nOpt * sum((A2MdlPrices(p)' - A2MktPrices).^2));

A2MdlImpv   = @(p) [A2MdlImpvCall(p); A2MdlImpvPut(p)];
% A2RMSE      = @(p) sqrt(1/A2nOpt * sum((A2MdlImpv(p) - A2mktvols).^2));

% MdlPrices       = @(p) [A1MdlPrices(p), A2MdlPrices(p)];
% 
% % Mkt Prices
% MktPrices       = [A1MktPrices; A2MktPrices];

% calibration procedure
qLB  = 0.2;
qUB  = 5;
options = optimoptions('fmincon', 'Display','off');

aux    = @(p) sqrt((1/(A1nOpt+A2nOpt)) * sum((MdlPrices(p)' - MktPrices).^2));
% qcal   = fmincon(@(p) aux(p), 1, [],[], [], [], qLB, qUB, [], options);
qcal   = fmincon(@(p) A1RMSE(p) + A2RMSE(p), 1, [],[], [], [], qLB, qUB, [], options);

%% Visualize the results
A1AvgMAPEvol = zeros(length(A1TTM),1);
A2AvgMAPEvol = zeros(length(A2TTM),1);
A1AvgMAPEprices = zeros(length(A1TTM),1);
A2AvgMAPEprices = zeros(length(A2TTM),1);


for i = 1:length(A1TTM)
    
    A1CalPricesCall = real(CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal1, Params, alpha, qcal));
    % A1CalPricesCall = real(CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1CallMoneyness{i}, A1TTM(i), cal1, Params, alpha));
    A1CalPricesPut  = real(CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal1, Params, alpha, qcal)) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})';
    % A1CalPricesPut  = real(CallPricesSatoNIGFFT(A1Forwards(i), A1Discounts(i), A1PutMoneyness{i}, A1TTM(i), cal1, Params, alpha)) - A1Discounts(i).*(A1Forwards(i) - A1PutStrikes{i})';
    A1mktvolsCall   = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CallPrices{i});
    A1CalVolCall    = blkimpv(A1Forwards(i), A1CallStrikes{i}, A1rates(i), A1TTM(i), A1CalPricesCall');
    A1mktvolsPut    = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1PutPrices{i}, 'Class','put');
    A1CalVolPut     = blkimpv(A1Forwards(i), A1PutStrikes{i}, A1rates(i), A1TTM(i), A1CalPricesPut', 'Class', 'put');
    A1MktPricesTemp = [A1CallPrices{i}; A1PutPrices{i}];
    A1mktvols       = [A1mktvolsCall; A1mktvolsPut];
    A1CalVol        = [A1CalVolCall; A1CalVolPut];
    A1CalPrices     = [A1CalPricesCall'; A1CalPricesPut'];
    Moneyness       = [-A1CallMoneyness{i}; -A1PutMoneyness{i}];
    A1AvgMAPEvol(i)    = sqrt(1/length(A1TTM).*sum(abs(A1mktvols - A1CalVol)./A1mktvols));
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
    ylim([0.1 0.3])
    titleName = strcat({'NIGSSD Smile for'}, {' '}, AssetNames{1,1}, {' '}, {'at'}, {' '}, {datestr(A1Maturities(i))});
    title(titleName, 'FontSize', FntSz, 'FontName', FntNm)
    xlabel('Moneyness', 'FontSize', FntSz, 'FontName', FntNm)
    ylabel('Implied Volatility', 'FontSize', FntSz, 'FontName', FntNm)
    grid on
    hold off

end

for i = 1:length(A2TTM)
    A2CalPricesCall = real(CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal2, Params, alpha, qcal));
    % A2CalPricesCall = real(CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2CallMoneyness{i}, A2TTM(i), cal2, Params, alpha));
    A2CalPricesPut  = real(CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal2, Params, alpha, qcal)) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})';
    % A2CalPricesPut  = real(CallPricesSatoNIGFFT(A2Forwards(i), A2Discounts(i), A2PutMoneyness{i}, A2TTM(i), cal2, Params, alpha)) - A2Discounts(i).*(A2Forwards(i) - A2PutStrikes{i})';
    A2mktvolsCall   = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CallPrices{i});
    A2CalVolCall    = blkimpv(A2Forwards(i), A2CallStrikes{i}, A2rates(i), A2TTM(i), A2CalPricesCall');
    A2mktvolsPut    = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2PutPrices{i}, 'Class','put');
    A2CalVolPut     = blkimpv(A2Forwards(i), A2PutStrikes{i}, A2rates(i), A2TTM(i), A2CalPricesPut', 'Class', 'put');
    A2MktPricesTemp = [A2CallPrices{i}; A2PutPrices{i}];
    A2mktvols       = [A2mktvolsCall; A2mktvolsPut];
    A2CalVol        = [A2CalVolCall; A2CalVolPut];
    A2CalPrices     = [A2CalPricesCall'; A2CalPricesPut'];
    A2AvgMAPEvol(i)    = sqrt(1/length(A2TTM).*sum(abs(A2mktvols - A2CalVol)./A2mktvols));
    A2AvgMAPEprices(i) = 100.*(1/length(A2MktPricesTemp).*sum(abs((A2MktPricesTemp-A2CalPrices))./A2MktPricesTemp));

    figure
    ax = gca;
    plot([-A2CallMoneyness{i}; -A2PutMoneyness{i}], A2mktvols, '+', 'MarkerSize', 10, 'DisplayName','MarketVols')
    hold on
    plot([-A2CallMoneyness{i}; -A2PutMoneyness{i}], A2CalVol, 'square' , 'MarkerSize', 10,'DisplayName','CalibratedVols')
    ax.XAxis.FontSize = FntSz;
    ax.YAxis.FontSize = FntSz;
    ax.XAxis.FontName = FntNm;
    ax.YAxis.FontName = FntNm;
    ylim([0.05 0.3])
    legend
    titleName = strcat({'NIGSSD Smile for '}, AssetNames{2,1}, {' '}, {'at'}, {' '}, {datestr(A2Maturities(i))});    
    title(titleName, 'FontSize', FntSz, 'FontName', FntNm)
    xlabel('Moneyness', 'FontSize', FntSz, 'FontName', FntNm)
    ylabel('Implied Volatility', 'FontSize', FntSz, 'FontName', FntNm)
    grid on
    hold off

end
A1CalPrices = A1MdlPrices(qcal)';
A2CalPrices = A2MdlPrices(qcal)';

%% Output
% MarginalParams = struct;
% MarginalParams.Asset1 = [A1a; A1beta; A1delta];
% MarginalParams.Asset2 = [A2a; A2beta; A2delta];
MarginalParams = [A1delta, A1theta, A1k; A2delta, A2theta, A2k];
CalibError = struct;
CalibError.volRMSE.Asset1 = sqrt(1/A1nOpt.*sum((A1CalVol - A1mktvols).^2));
CalibError.volRMSE.Asset2 = sqrt(1/A2nOpt.*sum((A2CalVol - A2mktvols).^2));
CalibError.PricesRMSE.Asset1 = sqrt(1/A1nOpt.*sum((A1CalPrices - A1MktPrices).^2));
CalibError.PricesRMSE.Asset2 = sqrt(1/A2nOpt.*sum((A2CalPrices - A2MktPrices).^2));
CalibError.volMAPE.Asset1 = A1AvgMAPEvol;
CalibError.volMAPE.Asset2 = A2AvgMAPEvol;
CalibError.PricesMAPE.Asset1 = A1AvgMAPEprices;
CalibError.PricesMAPE.Asset2 = A2AvgMAPEprices;

end