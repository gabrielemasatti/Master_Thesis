function [rho, a, RMSE, CalCorr, CorrBounds] = calibrateCorr(Asset1Params, Asset2Params, HistCorrelations, TimeHorizons, AssetNames, q)
%
% Function that calibrates the common parameters [rho, a] 
% for a pair of assets, matching the histCorrelations on different time
% horizons
%
% INPUT
% Asset1Params:     marginal params of first asset  [a, beta, delta]
% Asset2Params:     marginal params of second asset
% HistCorrelations: historical corr on TimeHorizons
% TimeHorizons:     time horizons where the hist correlations are computed
% AssetNames:       names of the couple of assets considered
% q:                Sato exponent considered
%
% OUTPUT
% rho:              rho_{ij}
% a:                a common parameter
% RMSE:             average RMSE on hist correlations
% CalCorr:          Calibrated values of the correlations
%

% Preliminaries
FntSz = 20;
FntNm = 'Times';
TimeHorizons = TimeHorizons./365; % fraction of year (approx)

%% Calibration of a, rho on log-returns with different windows
options = optimoptions("fmincon", "OptimalityTolerance", 1e-10);
x00 = [0.2; 0.2];           % [rho, a]
aux1 = @(rho, a) sqrt((1/size(TimeHorizons, 1)).*sum((ComputeCorrelation(Asset1Params, Asset2Params, TimeHorizons, rho, a, q) - HistCorrelations).^2));
CommonParams = fmincon(@(cal) aux1(cal(1), cal(2)), x00, [], [], [], [], [-1; 0], [1-1e-2; min(1/(sqrt(Asset1Params(1))), 1/(sqrt(Asset2Params(1)))) - 1e-2], [], options);
% CommonParams = lsqnonlin(@(cal) arrayfun(@(i) (ComputeCorrelation(Asset1Params, Asset2Params, TimeHorizons(i), cal(1), cal(2), q) - HistCorrelations(i)).^2, [1:length(TimeHorizons)]'), x00, [-1; 0], [1; min(1/(sqrt(Asset1Params(1))), 1/(sqrt(Asset2Params(1))))]);

rho = CommonParams(1);
a   = CommonParams(2);

%% Boundaries of correlations

Zmean  = sqrt(pi)*a;
ZVar   = 0.5*sqrt(pi)*a;
S1mean = sqrt(pi*Asset1Params(1));
S1Var  = 0.5*sqrt(pi*Asset1Params(1)^3);
S2mean = sqrt(pi*Asset2Params(1));
S2Var  = 0.5*sqrt(pi*Asset2Params(1)^3);

rhoMIN = rho*sqrt(Asset1Params(1)*Asset2Params(1))*Zmean/sqrt(S1mean*S2mean);
rhoMAX = (Asset1Params(2)*Asset2Params(2)*Asset1Params(1)*Asset2Params(1)*ZVar)/(abs(Asset1Params(2))*abs(Asset2Params(2))*sqrt(S1Var*S2Var));

CorrBounds = [rhoMIN; rhoMAX];

% plot of the historical corr wrt bounds
% figure()
% plot(TimeHorizons, rhoMIN*ones(size(TimeHorizons,1),1), 'LineWidth', 2,'DisplayName', 'LB')
% hold on
% plot(TimeHorizons, rhoMAX*ones(size(TimeHorizons,1),1), 'LineWidth', 2, 'DisplayName', 'UB')
% hold on
% plot(TimeHorizons, HistCorrelations, '-*', 'LineWidth', 2, 'DisplayName', 'HistCorrelation')
% legend
% title(['Correlation Bounds vs HistCorr', ' ', AssetNames{1,1}, ' ', AssetNames{2,1}], FontName=FntNm, FontSize=FntSz)
% xlabel('TimeHorizons', FontName=FntNm, FontSize=FntSz)
% grid on

% output values
CalCorr = ComputeCorrelation(Asset1Params, Asset2Params, TimeHorizons, rho, a, q); 
RMSE    = sqrt((1/size(TimeHorizons,1)).*sum((ComputeCorrelation(Asset1Params, Asset2Params, TimeHorizons, rho, a, q) - HistCorrelations).^2));

end