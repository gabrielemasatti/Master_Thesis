function [Correlations] = AssetCorrs(Prices, TimeHorizons, AssetNames)
% 
% Returns the correlations of the assets with different horizons computed
% on different windows
% 
% INPUT
% Prices:        matrix containig the series of prices with dates
% AssetNames:    Names of the assets of which the corr is wanted
% TimeHorizons:  Time horizons on which correlation is wanted
% 
% OUTPUT
% Correlations: matrix with [Time horizons, Correlation of the log returns]
%
% CALLS
% Corr
%

FntNm = 'Times';
FntSz = 20;

% Preallocate the matrix
FwdCorrelations = zeros(size(TimeHorizons,1), 1);
BwdCorrelations = zeros(size(TimeHorizons,1), 1);


% Directly load from excel
PriceDates = load("PriceDates.mat");
% PriceDates = load("SPXN225PriceDates.mat");
PriceDates = PriceDates.Dates;

% Actual computation of returns and correlations

for i=1:size(TimeHorizons)
    
    FwdSPXReturns       = log(Prices(TimeHorizons(i)+1:TimeHorizons(i):end, 1)./Prices(1:TimeHorizons(i):end-TimeHorizons(i), 1));
    FwdSX5EReturns      = log(Prices(TimeHorizons(i)+1:TimeHorizons(i):end, 2)./Prices(1:TimeHorizons(i):end-TimeHorizons(i), 2));
    FwdCorrelations(i)  = corr(FwdSPXReturns, FwdSX5EReturns);

    BwdSPXReturns       = log(Prices(end:-TimeHorizons(i):1+TimeHorizons(i), 1)./Prices(end-TimeHorizons(i):-TimeHorizons(i):1, 1));
    BwdSX5EReturns      = log(Prices(end:-TimeHorizons(i):1+TimeHorizons(i), 2)./Prices(end-TimeHorizons(i):-TimeHorizons(i):1, 2));
    BwdCorrelations(i)  = corr(BwdSPXReturns, BwdSX5EReturns);
    
    
%     figure
%     ax = gca;
%     plot(PriceDates(TimeHorizons(i)+1:TimeHorizons(i):end), SX5EReturns, 'DisplayName', AssetNames{1,1})
%     hold on
%     plot(PriceDates(TimeHorizons(i)+1:TimeHorizons(i):end), SPXReturns, 'DisplayName', AssetNames{1,2})
%     datetick('x', 20)
%     legend
%     hold off
%     grid on
%     ax.XAxis.FontSize = 16;
%     ax.YAxis.FontSize = 16;
%     ax.XAxis.FontName = FntNm;
%     ax.YAxis.FontName = FntNm;

end

hold off

Correlations = [TimeHorizons, BwdCorrelations];

end