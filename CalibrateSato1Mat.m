function [MarginalParams, CommonParams, Errors, CalibCorr, CorrBounds, qCal] = CalibrateSato1Mat(RequiredMat, SettlementDate, ...
    A1Discounts, A2Discounts, A1Forwards, A2Forwards, A1Data, A2Data, q, alpha, x00, AssetNames, Corr)
%
% Calibrate the model in (2.1 of 11072023) 
% /considers two assets and calibrates their marginal params [a, beta,
% delta] and correlations on one single maturity
% 
% INPUT
% RequiredMat:    Date at which the marginals are calibrated
% SettlementDate: Date at which data are observed
% A1Discounts:    Discount factors for first asset
% A2Discounts:    Forward Prices for first asset
% A1Forwards:     Discount factors for second asset
% A2Forwards:     Forward Prices for second asset
% A1Data:         Option data for first asset
% A2Data:         Option data for second asset
% q:              Sato exponent
% alpha:          Parameter of EtaS
% x00:            Starting point for the calibration
% AssetNames:     Names of the financial assets considered
%
% OUTPUT
% MarginalParams: [a_j, beta_j, delta_j] for every marginal
% CommonParams:   [a, rho]
% Errors:         Struct containing calib errors 
% CalibCorr:      Calibrated correlation values for 1d, 5d, 20d, 120d
%
% CALLS
% filterOTMoptions
% CalibrateSatoToVolatilitySurfaceCallPut
% calibrateCorr
%

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-

%% Marginal Parameters
MarginalParams = zeros(2,3);

% Filter Options (only OTM) at required date
[A1OTMCall, A1OTMPut] = filterOTMoptions(RequiredMat, SettlementDate, A1Data, A1Forwards, A1Discounts, unique(A1Data.MATURITIES));
[A2OTMCall, A2OTMPut] = filterOTMoptions(RequiredMat, SettlementDate, A2Data, A2Forwards, A2Discounts, unique(A2Data.MATURITIES));

% Find right indices
A1IdxDate = find(~(A1OTMCall.MATURITIES(1) - unique(A1Data.MATURITIES)));
A1nCall   = size(A1OTMCall,1);
A2IdxDate = find(~(A2OTMCall.MATURITIES(1) - unique(A2Data.MATURITIES)));
A2nCall   = size(A2OTMCall,1);

% Calibration on the specified date
[A1a, A1beta, A1delta, A1RMSE, A1MAPE] = CalibrateSatoToVolatilitySurfaceCallPut(SettlementDate, ...
    A1Discounts(A1IdxDate), A1Forwards(A1IdxDate), [A1OTMCall.STRIKES; A1OTMPut.STRIKES], A1nCall,...
    [A1OTMCall.MID; A1OTMPut.MID], A1OTMCall.MATURITIES(1), x00, q, alpha, AssetNames{1,1});

[A2a, A2beta, A2delta, A2RMSE, A2MAPE] = CalibrateSatoToVolatilitySurfaceCallPut(SettlementDate, ...
    A2Discounts(A2IdxDate), A2Forwards(A2IdxDate), [A2OTMCall.STRIKES; A2OTMPut.STRIKES], A2nCall,...
    [A2OTMCall.MID; A2OTMPut.MID], A2OTMCall.MATURITIES(1), x00,  q, alpha, AssetNames{1,2});

% Structs with info on both the surfaces
Asset1Strct                = struct;
Asset1Strct.Discount       = A1Discounts(A1IdxDate);
Asset1Strct.Forward        = A1Forwards(A1IdxDate);
Asset1Strct.CallStrikes    = A1OTMCall.STRIKES;
Asset1Strct.PutStrikes     = A1OTMPut.STRIKES;
Asset1Strct.nCall          = A1nCall;
Asset1Strct.CallMid        = A1OTMCall.MID;
Asset1Strct.PutMid         = A1OTMPut.MID;
Asset1Strct.Maturity       = A1OTMCall.MATURITIES(1);
Asset1Strct.MarginalParams = [A1a; A1beta; A1delta];

Asset2Strct                = struct;
Asset2Strct.Discount       = A2Discounts(A2IdxDate);
Asset2Strct.Forward        = A2Forwards(A2IdxDate);
Asset2Strct.CallStrikes    = A2OTMCall.STRIKES;
Asset2Strct.PutStrikes     = A2OTMPut.STRIKES;
Asset2Strct.nCall          = A2nCall;
Asset2Strct.CallMid        = A2OTMCall.MID;
Asset2Strct.PutMid         = A2OTMPut.MID;
Asset2Strct.Maturity       = A2OTMCall.MATURITIES(1);
Asset2Strct.MarginalParams = [A2a; A2beta; A2delta];

% Calibration of q
qCal = CalibrateQParam(SettlementDate, Asset1Strct, Asset2Strct, alpha, AssetNames);

% output
MarginalParams(1,:) = [A1a, A1beta, A1delta];
MarginalParams(2,:) = [A2a, A2beta, A2delta];
Errors.Marginals.RMSE(1) = A1RMSE;
Errors.Marginals.RMSE(2) = A2RMSE;
Errors.Marginals.MAPE(1) = A1MAPE;
Errors.Marginals.MAPE(2) = A2MAPE;

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-

%% Common Parameters
[rho, a, RMSECorr, CalibCorr, CorrBounds] = calibrateCorr(MarginalParams(1,:), MarginalParams(2,:), Corr(:,2), Corr(:, 1), AssetNames, qCal);
CommonParams = [rho; a];
Errors.Common.RMSE = RMSECorr;

end