function [Price, Correlations] = MBRCWOFCOPULApricing(Nsim, ULsMktData, SetDate, Ncoupons, CouponFreq, CouponBarrier, CouponLevel, FinalBarrier, MarginalParams, rho, q)
%
% Function that provides the price of a MBRC with a WOF feature
%
% INPUT
% ULsMktData:     structure with market data of both assets
% SetDate:        settlement date
% Ncoupons:       n of coupons payable by the contract
% CouponFreq:     freq in months of the coupons
% MarginalParams: parameters of marginal distributions of the ULs
%
% OUTPUT
% Price: Price of the MBRCWOF obtained via montecarlo simulations
%
%
rng(0)

%% ULs parameters
A1S0        = ULsMktData.Asset1.S0;
A1Dates     = ULsMktData.Asset1.Dates;
A1Forwards  = ULsMktData.Asset1.Forwards;
A1Discounts = ULsMktData.Asset1.Discounts;
A1Rates     = ULsMktData.Asset1.Rates;
A1div       = ULsMktData.Asset1.Div;

A2S0        = ULsMktData.Asset2.S0;
A2Dates     = ULsMktData.Asset2.Dates;
A2Forwards  = ULsMktData.Asset2.Forwards;
A2Discounts = ULsMktData.Asset2.Discounts;
A2Rates     = ULsMktData.Asset2.Rates;
A2div       = ULsMktData.Asset2.Div;

%% ULs Rates, dividends and Product Dates
MonDates    = zeros(Ncoupons, 1);    % last is the maturity (first is setdate) 

for i = 1:Ncoupons    
    MonDates(i)   = busdate(datenum(addtodate(SetDate, CouponFreq*(i), 'month')));
end
MonDates = [SetDate; MonDates];


% interpolation of rates
A1MonRates = interp1(A1Dates, A1Rates, MonDates(2:end));
A1MonDiv   = interp1(A1Dates, A1div, MonDates(2:end));
A1MonDisc  = exp(-A1MonRates.*yearfrac(SetDate, MonDates(2:end), 3));
A1MonDiscDiv    = exp(-A1MonDiv.*yearfrac(SetDate, MonDates(2:end), 3));
A1FwdMonDisc    = A1MonDisc(2:end)./A1MonDisc(1:end-1);
A1FwdMonDiscDiv = A1MonDiscDiv(2:end)./A1MonDiscDiv(1:end-1);
A1FwdMonRates   = - log(A1FwdMonDisc)./yearfrac(MonDates(2:end-1), MonDates(3:end), 3);
A1FwdMonDiv     = - log(A1FwdMonDiscDiv)./yearfrac(MonDates(2:end-1), MonDates(3:end), 3);

A2MonRates = interp1(A2Dates, A2Rates, MonDates(2:end));
A2MonDiv   = interp1(A2Dates, A2div, MonDates(2:end));
A2MonDisc  = exp(-A2MonRates.*yearfrac(SetDate, MonDates(2:end), 3));
A2MonDiscDiv    = exp(-A2MonDiv.*yearfrac(SetDate, MonDates(2:end), 3));
A2FwdMonDisc    = A2MonDisc(2:end)./A2MonDisc(1:end-1);
A2FwdMonDiscDiv = A2MonDiscDiv(2:end)./A2MonDiscDiv(1:end-1);
A2FwdMonRates   = - log(A2FwdMonDisc)./yearfrac(MonDates(2:end-1), MonDates(3:end), 3);
A2FwdMonDiv     = - log(A2FwdMonDiscDiv)./yearfrac(MonDates(2:end-1), MonDates(3:end), 3);

% interpolation @ maturity
% A1rateMat  = interp1(A1Dates, A1Rates, MatDate);
% A1DiscMat  = exp(-A1rateMat.*yearfrac(SetDate, MatDate));
% 
% A2rateMat  = interp1(A2Dates, A2Rates, MatDate);
% A2DiscMat  = exp(-A2rateMat.*yearfrac(SetDate, MatDate));
% 
% A1divMat   = interp1(A1Dates, A1div, MatDate);
% A2divMat   = interp1(A2Dates, A2div, MatDate);

% initial forward @ first mon date
A1Fwd = A1S0*exp((A1MonRates(1) - A1MonDiv(1))*yearfrac(SetDate, MonDates(2), 3));
A2Fwd = A2S0*exp((A2MonRates(1)- A1MonDiv(1))*yearfrac(SetDate, MonDates(2), 3));

% initial forward @ mat
% A1Fwd = A1S0*exp(A1rateMat - A1divMat)*TTM;
% A2Fwd = A2S0*exp(A2rateMat - A2divMat)*TTM;

% preallocate variables
f1 = zeros(Nsim, Ncoupons);
F1path = [A1Fwd*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
S1path = [A1S0*ones(Nsim,1), zeros(Nsim, Ncoupons)];

f2 = zeros(Nsim, Ncoupons);
F2path = [A2Fwd*ones(Nsim,1), zeros(Nsim, Ncoupons-1)];
S2path = [A2S0*ones(Nsim,1), zeros(Nsim, Ncoupons)];

Perf1  = zeros(Nsim, Ncoupons);
Perf2  = zeros(Nsim, Ncoupons);

WOFPerf = zeros(Nsim, Ncoupons);
DigPayout  = zeros(Nsim, Ncoupons);
CouponProb   = zeros(Ncoupons, 1);
DigDiscPayout = zeros(Nsim, Ncoupons);
Correlations = [corr(A1S0*ones(Nsim,1), A2S0*ones(Nsim,1)); zeros(Ncoupons, 1)];

A1LogReturnsCompounded = zeros(Nsim, Ncoupons);
A2LogReturnsCompounded = zeros(Nsim, Ncoupons);

%% simulations and Digital Payout 
for i = 1:Ncoupons

    t = yearfrac(SetDate, MonDates(i+1), 3);
    s = yearfrac(SetDate, MonDates(i), 3);
    
    [f1(:,i), f2(:,i)] = simulateNIGCOPULA(Nsim, t, s, MarginalParams, rho, q);
    S1path(:, i+1) = F1path(:, i).*exp(f1(:,i));     % F1(t_i,t_i) = S1(t_i)
    S2path(:, i+1) = F2path(:, i).*exp(f2(:,i));     % F2(t_i,t_i) = S2(t_i)

    Perf1(:, i) = S1path(:,i+1)./S1path(:,1);
    Perf2(:, i) = S2path(:,i+1)./S2path(:,1);
    
    if (i~=Ncoupons)
        F1path(:, i+1) = S1path(:, i+1).*exp((A1FwdMonRates(i)-A1FwdMonDiv(i))*yearfrac(MonDates(i), MonDates(i+1), 3));
        F2path(:, i+1) = S2path(:, i+1).*exp((A2FwdMonRates(i)-A2FwdMonDiv(i))*yearfrac(MonDates(i), MonDates(i+1), 3));
    end

    WOFPerf(:, i)      = min(Perf1(:, i), Perf2(:, i));

    DigPayout(:, i) = (CouponLevel.*ones(Nsim, 1)).*(WOFPerf(:,i) >= CouponBarrier);

    if (i~=1)
        DigPayout(:, i) = DigPayout(:, i).*(DigPayout(:, i-1) > 0);
    end

    DigDiscPayout(:, i) = A1MonDisc(i).*DigPayout(:, i);

%     if (i~=Ncoupons)
%         DiscPayout(:, i) = A1MonDisc(i).*Payout(:, i);
%     end
    
    CouponProb(i) = sum(DigPayout(:, i) > 0)./Nsim;

end

Correlations(end) = corr(log(S1path(:, end)./S1path(:, 1)), log(S2path(:, end)./S2path(:, 1)));

% Correlation check
TempDiff  = (yearfrac(SetDate, MonDates(2:end), 3).^q - yearfrac(SetDate, MonDates(1:end-1), 3).^q).^2./(yearfrac(SetDate, MonDates(2:end), 3).^(2*q));
TempRatio = (yearfrac(SetDate, MonDates(2:end-1), 3)./yearfrac(SetDate, MonDates(3:end), 3)).^(2*q);
TotProd   = cumprod(TempRatio);
TotSum    = TempDiff(1)*TotProd(end);
rho       = 0.7475;

for i = 2:length(TempDiff)-1
    TotProd = cumprod(TempRatio(i:end));
    TotSum = TotSum + TempDiff(i)*TotProd(end);
end

TotCorr = rho*(((yearfrac(SetDate, MonDates(end), 3).^q - yearfrac(SetDate, MonDates(end-1), 3).^q)).^2/(yearfrac(SetDate, MonDates(end), 3).^(2*q)) + ...
   TotSum);

% Payout at maturity

MatPayout = ones(Nsim, 1).*(WOFPerf(:, end) >= FinalBarrier)+ ...
    WOFPerf(:, end).*(WOFPerf(:, end) < FinalBarrier);

MatDiscPayout = A1MonDisc(end).*MatPayout;

DiscPayout = MatDiscPayout + sum(DigDiscPayout,2);

% Price
DigPrice   = mean(sum(DigDiscPayout,2));
ZCBDIPrice = mean(MatDiscPayout);

% Price = mean(DiscPayout);
Price = DigPrice + ZCBDIPrice;














end