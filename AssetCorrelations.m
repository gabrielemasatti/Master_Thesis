function [Correlations, Returns] = AssetCorrelations(DataFile)
% 
% Returns the correlations of the assets with different horizons computed
% on different windows
% 
% INPUT
% DataFile:     Name of the file where the prices of the assets are stored
% NB: the dataset must contain already filtered prices, with non overlapping 
% periods for longer than 1 day horizons, stored in different sheets
%
% OUTPUT
% Correlations: matrix with Correlation of the log returns
% Returns:      struct with daily, weekly, monthly, quarterly and
%               semiannually log returns
% 
% CALLS
% Corr
%

% Preliminaries
FntNm = 'Times';
FntSz = 20;

% Store the vector of prices
DailyPrices        = xlsread(DataFile, "Daily");
WeeklyPrices       = xlsread(DataFile, "Weekly");
MonthlyPrices      = xlsread(DataFile, "Monthly");
QuarterlyPrices    = xlsread(DataFile, "Quarterly");
SemiAnnuallyPrices = xlsread(DataFile, "SemiAnnually");
AnnuallyPrices     = xlsread(DataFile, "Annually");

% Returns and Historical Correlations
A1DailyReturns         = log(DailyPrices(2:end, 2)./DailyPrices(1:end-1, 2));
A2DailyReturns         = log(DailyPrices(2:end, 3)./DailyPrices(1:end-1, 3));
DailyCorr              = corr(A1DailyReturns, A2DailyReturns);

A1WeeklyReturns        = log(WeeklyPrices(2:end, 2)./WeeklyPrices(1:end-1, 2));
A2WeeklyReturns        = log(WeeklyPrices(2:end, 3)./WeeklyPrices(1:end-1, 3));
WeeklyCorr             = corr(A1WeeklyReturns, A2WeeklyReturns);

A1MonthlyReturns       = log(MonthlyPrices(2:end, 2)./MonthlyPrices(1:end-1, 2));
A2MonthlyReturns       = log(MonthlyPrices(2:end, 3)./MonthlyPrices(1:end-1, 3));
MonthlyCorr            = corr(A1MonthlyReturns, A2MonthlyReturns);

A1QuarterlyReturns     = log(QuarterlyPrices(2:end, 2)./QuarterlyPrices(1:end-1, 2));
A2QuarterlyReturns     = log(QuarterlyPrices(2:end, 3)./QuarterlyPrices(1:end-1, 3));
QuarterlyCorr          = corr(A1QuarterlyReturns, A2QuarterlyReturns);

A1SemiAnnuallyReturns = log(SemiAnnuallyPrices(2:end, 2)./SemiAnnuallyPrices(1:end-1, 2));
A2SemiAnnuallyReturns = log(SemiAnnuallyPrices(2:end, 3)./SemiAnnuallyPrices(1:end-1, 3));
SemiAnnuallyCorr      = corr(A1SemiAnnuallyReturns, A2SemiAnnuallyReturns);

A1AnnuallyReturns = log(AnnuallyPrices(2:end, 2)./AnnuallyPrices(1:end-1, 2));
A2AnnuallyReturns = log(AnnuallyPrices(2:end, 3)./AnnuallyPrices(1:end-1, 3));
AnnuallyCorr      = corr(A1AnnuallyReturns, A2AnnuallyReturns);


% Correlations for the output
% Correlations = [DailyCorr; WeeklyCorr; MonthlyCorr; QuarterlyCorr; SemiAnnuallyCorr];
Correlations = [DailyCorr; WeeklyCorr; MonthlyCorr; QuarterlyCorr; SemiAnnuallyCorr; AnnuallyCorr];

% Returns Struct
Returns              = struct;
Returns.Daily        = [A1DailyReturns, A2DailyReturns];
Returns.Weekly       = [A1WeeklyReturns, A2WeeklyReturns];
Returns.Monthly      = [A1MonthlyReturns, A2MonthlyReturns];
Returns.Quarterly    = [A1QuarterlyReturns, A2QuarterlyReturns];
Returns.SemiAnnually = [A1SemiAnnuallyReturns, A1SemiAnnuallyReturns];
Returns.Annually     = [A1AnnuallyReturns, A2AnnuallyReturns];

end