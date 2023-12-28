function [a, beta, delta, RMSE, MAPE] = CalibrateSatoToVolatilitySurfaceCallPut(setDate, discount, Forward, strikes, nCall, prices, maturity, x00, q, alpha, AssetName)
%
% Function that calibrates the parameters of a Sato NIG
% process at a certain maturity, making use of OTM calls & puts whose prices are
% stored in prices. 
% It performs also a plot of implied volatility smile obtained and computes the most common errors
% MSE & MAPE
%  
% INPUT
% setDate:   settlement date considered
% discount:  discount factor @maturity
% Forward:   forward price with expiry @maturity (F(t0,t))
% strikes:   strikes of the options
% nCall:     number of calls between options provided
% prices:    prices of options passed as input
% maturity:  maturity of the options
% q:         parameter of ATS process (0.5 -> NIG, 0 -> VG)
% x00:       starting point of the calibration
%
% OUTPUT
% a:       calibrated model parameter which represents the avarage volatility
% beta:    calibrated model parameter which represents 'volatility of the volatility'
% delta:   calibrated model parameter which represents the skewness of the volatility curve
% MSE:     Mean Squared Error of the performed calibration
% MAPE:    Mean Absolute Percentage Error of the performed calibration
% mktSkew: implied skew of the market (difference between last & first impl vols from market prices)
% calSkew: reproduced implied skew  (difference between last & first impl vols from calibrated prices)
%
% CALLS 
% fmincon    for the constrained minimization of the distance (NIG constraint only)
% CallPricesNMVMFFT
% 

%% Displaying data
IBdaycount = 3;
TTM        = yearfrac(setDate, maturity, IBdaycount);
moneyness  = log(Forward./strikes);
rate       = -log(discount)/TTM;

%% Model Prices
% Calibrating  parameters for the FFT
I_res     = 2*pi*exp(-sign(moneyness)*0.5.*moneyness);
M         = 15;
options   = optimset('TolFun',1e-5);
x0        = 0.0025;
LB        = eps;
UB        = 0.01;
fTS       = @(v) 1./(v.^2 + 0.25);
ff        = @(x) 1./(x.^2 + 1/4);
dz        = lsqnonlin(@(dz) abs( FourierTransform(fTS, moneyness,  M, dz)- I_res), x0, LB, UB, options); % find dz such that FFT replicates the residual integral
Params    = FFTparameters(M, dz, 1);
I         = computeIntegral(ff, moneyness, [], Params, 1);
MdlPricesCall = @(cal) CallPricesSatoFFT(Forward, discount, moneyness(1:nCall), TTM, cal, Params, alpha, q);
MdlPricesPut  = @(cal) CallPricesSatoFFT(Forward, discount, moneyness(nCall+1:end), TTM, cal, Params, alpha, q) - discount.*(Forward - strikes(nCall+1:end))'; 
MdlPrices     = @(cal) [MdlPricesCall(cal), MdlPricesPut(cal)];


% 
% errorFFTINT = abs(I - I_res);
% figure()
% plot(moneyness, I_res, '*b', 'LineWidth', 2)
% hold on
% plot(moneyness, I, '+g', 'LineWidth', 2)
% grid on
% text(moneyness - 0.0002, min(I, I_res) , num2str(abs(I - I_res)),'Color', 'g', 'FontSize', 8); % error for each point 
% tINT = text(min(moneyness), I(2), ['\bf Error : ', num2str(sum(errorFFTINT))], 'Color', 'g');  % total error
% tINT.FontSize = 13;
% legend('Residuals', 'FFT')
% hold off

% NIG constraints (passed to fmincon)
% Matlab help : easy way to implement a non linear constraint 
% function [c,ceq] = constraints (cal)
%   w   = 1/(2*cal(3).*cal(1).^2);  
%   c   = [-w-cal(2),...
%       -cal(1)];
%   ceq = [];
% end

% starting points, lower bounds & upper bounds
% x00 = [0.3; 3; TTM];     % sigma, eta, k
% LB = [0.17; 8; 0.2];
LB = [0; -30; 0.0001];   % a, beta, delta
% UB = [1; 20; 5];
UB = [100; 60; 0.5];

aux = @(cal) 1/length(prices)*sqrt(sum(abs(prices-MdlPrices(cal)').^2)); % function that we want to minimize the weigths are negligible since they are equal to 1

% implied volatility calibration
% aux             = @(cal) sqrt(sum(abs(arrayfun(@(i) blkimpv(Forward, strikes(i), rate, TTM, prices(i)) ...
%     -blkimpv(Forward, strikes(i), rate, TTM, MdlPrices(cal))').^2));

%% Calibration & Results
% const = @constraints;

% minimization
% cal   = fmincon(@(cal) aux(cal),x00,[],[],[],[],LB,UB,const);
cal   = fmincon(@(cal) aux(cal),x00,[],[],[],[],LB,UB);

% calibrated parameters
a     = cal(1);
beta  = cal(2);
delta = cal(3);
gamma = sqrt((1/a+delta^2*beta^2)*1/(delta^2));

%% Compute Model Prices with calibrated parameters
CalPricesCall = real(CallPricesSatoFFT(Forward, discount, moneyness(1:nCall), TTM, cal, Params, alpha, q));
CalPricesPut  = real(CallPricesSatoFFT(Forward, discount, moneyness(nCall+1:end), TTM, cal, Params, alpha, q)) - discount.*(Forward - strikes(nCall+1:end))';
CalPrices     = [CalPricesCall, CalPricesPut];

%% Black Volatility obtain from the Calibrated Model Prices 
mktvolsCall   = blkimpv(Forward, strikes(1:nCall), rate, TTM, prices(1:nCall));
CalVolCall    = blkimpv(Forward, strikes(1:nCall), rate, TTM, CalPricesCall');
mktvolsPut    = blkimpv(Forward, strikes(nCall+1:end), rate, TTM, prices(nCall+1:end), 'Class','put');
CalVolPut     = blkimpv(Forward, strikes(nCall+1:end), rate, TTM, CalPricesPut', 'Class', 'put');
mktvols       = [mktvolsCall; mktvolsPut];
CalVol        = [CalVolCall; CalVolPut];

% calibration errors (on prices & volatilities)
PriceSquaredErrors = sqrt(abs(prices - CalPrices').^2);
PercError          = abs(prices - CalPrices')./prices;
MAPE               = 100*mean(PercError);
RMSE               = mean(PriceSquaredErrors);
moneyness          = log(strikes./Forward);

% plot
FntNm = 'Times';
FntSz = 20;
figure()
ax = gca;
plot(moneyness, CalVol,'+', 'MarkerSize', 8)
hold on 
plot(moneyness, mktvols, 'square', 'MarkerSize',8)
grid on
legend('Model Volatilities', 'Market Volatilities')
xlabel('Moneyness', FontName=FntNm, FontSize=FntSz)
ylabel('Volatilities', FontName=FntNm, FontSize=FntSz)
title([AssetName, ' ', 'Smile S-IG Sato @', num2str(datestr(maturity))], FontName=FntNm, FontSize=FntSz)
ax.XAxis.FontSize = FntSz;
ax.YAxis.FontSize = FntSz;
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;

end