function qCal = CalibrateQParam(SettlementDate, A1Strct, A2Strct, alpha, AssetNames)
%
% Function that given the calibrated marginal parameters on one maturity,
% calibrates the common q, minimizing the global RMSE
% 
% INPUT
% SettlementDate: SettlementDate
% A1Struct:       struct containing all relevant info on Asset1
% A2Struct:       struct containing all relevant info on Asset2
% x00:            starting point of q
% alpha:          alpha parameter of the ETaS
% AssetNames:     names of the assets
% 
% OUTPUT
% qCal:           Calibrated common q parameter
%

%% Displaying data for Asset1
IBdaycount      = 3;
A1TTM           = yearfrac(SettlementDate, A1Strct.Maturity, IBdaycount);
A1Strct.Strikes = [A1Strct.CallStrikes; A1Strct.PutStrikes];
A1moneyness     = log(A1Strct.Forward./A1Strct.Strikes);
A1rate          = -log(A1Strct.Discount)/A1TTM;
A1CallMktPrices = A1Strct.CallMid;
A1PutMktPrices  = A1Strct.PutMid;
A1Prices        = [A1CallMktPrices; A1PutMktPrices];


%% Model Prices for Asset1
% Calibrating  parameters for the FFT
I_res     = 2*pi*exp(-sign(A1moneyness)*0.5.*A1moneyness);
M         = 15;
options   = optimset('TolFun',1e-5);
x0        = 0.0025;
LB        = eps;
UB        = 0.01;
fTS       = @(v) 1./(v.^2 + 0.25);
ff        = @(x) 1./(x.^2 + 1/4);
dz        = lsqnonlin(@(dz) abs(FourierTransform(fTS, A1moneyness,  M, dz)- I_res), x0, LB, UB, options); % find dz such that FFT replicates the residual integral
Params    = FFTparameters(M, dz, 1);
I         = computeIntegral(ff, A1moneyness, [], Params, 1);
A1MdlPricesCall = @(q) CallPricesSatoFFT(A1Strct.Forward, A1Strct.Discount, A1moneyness(1:A1Strct.nCall), A1TTM, A1Strct.MarginalParams, Params, alpha, q);
A1MdlPricesPut  = @(q) CallPricesSatoFFT(A1Strct.Forward, A1Strct.Discount, A1moneyness(A1Strct.nCall+1:end), A1TTM, A1Strct.MarginalParams, Params, alpha, q) - A1Strct.Discount.*(A1Strct.Forward - A1Strct.Strikes(A1Strct.nCall+1:end))'; 
A1MdlPrices     = @(q) [A1MdlPricesCall(q), A1MdlPricesPut(q)];

%% Displaying data for Asset2
IBdaycount  = 3;
A2TTM       = yearfrac(SettlementDate, A2Strct.Maturity, IBdaycount);
A2Strct.Strikes = [A2Strct.CallStrikes; A2Strct.PutStrikes];
A2moneyness = log(A2Strct.Forward./A2Strct.Strikes);
A2rate      = -log(A2Strct.Discount)/A2TTM;
A2CallMktPrices = A2Strct.CallMid;
A2PutMktPrices  = A2Strct.PutMid;
A2Prices        = [A2CallMktPrices; A2PutMktPrices];

%% Model Prices for Asset2
% Calibrating  parameters for the FFT
I_res     = 2*pi*exp(-sign(A2moneyness)*0.5.*A2moneyness);
M         = 15;
options   = optimset('TolFun',1e-5);
x0        = 0.0025;
LB        = eps;
UB        = 0.01;
fTS       = @(v) 1./(v.^2 + 0.25);
ff        = @(x) 1./(x.^2 + 1/4);
dz        = lsqnonlin(@(dz) abs(FourierTransform(fTS, A2moneyness,  M, dz)- I_res), x0, LB, UB, options); % find dz such that FFT replicates the residual integral
Params    = FFTparameters(M, dz, 1);
I         = computeIntegral(ff, A2moneyness, [], Params, 1);
A2MdlPricesCall = @(q) CallPricesSatoFFT(A2Strct.Forward, A2Strct.Discount, A2moneyness(1:A2Strct.nCall), A2TTM, A2Strct.MarginalParams, Params, alpha, q);
A2MdlPricesPut  = @(q) CallPricesSatoFFT(A2Strct.Forward, A2Strct.Discount, A2moneyness(A2Strct.nCall+1:end), A2TTM, A2Strct.MarginalParams, Params, alpha, q) - A2Strct.Discount.*(A2Strct.Forward - A2Strct.Strikes(A2Strct.nCall+1:end))'; 
A2MdlPrices     = @(q) [A2MdlPricesCall(q), A2MdlPricesPut(q)];

% starting points, lower bounds & upper bounds
LB = 0.4;   
UB = 3;

% Minimize RMSE of both the surfaces
Prices    = [A1Prices; A2Prices];
MdlPrices = @(q) [A1MdlPrices(q), A2MdlPrices(q)];

aux  = @(q) sqrt(1/length(Prices)*sum(abs(Prices-MdlPrices(q)').^2)); % function that we want to minimize the weigths are negligible since they are equal to 1
qCal = fmincon(@(q) aux(q), 1, [], [], [], [], LB, UB);

end
