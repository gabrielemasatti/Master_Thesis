function [CallOneDateTab, PutOneDateTab] = filterOTMoptions(maturity, settlement, OptTab, Forwards, Discounts, Dates)
%
% Filters the options according to maturity and liquidity
%
% INPUT
% maturity
% settlement
% OptTab
% Forwards
% Discounts
% Dates
%
% OUTPUT
% 
%
%

IBDaycount    = 3;
TTMvec        = yearfrac(settlement, Dates, IBDaycount);
[~, matIndex] = min(abs(maturity - TTMvec));                  % indexes of the correct maturity
OptMat        = Dates(matIndex);                              % maturity as a datenum
TTM           = yearfrac(settlement, OptMat, IBDaycount);

OptIdx  = find(OptTab.MATURITIES == OptMat);

FiltOptTab  = OptTab(OptIdx, :);

Fwd      = Forwards(matIndex);
ZeroRate = -log(Discounts(matIndex))/TTM;

% consider only OTM Options for the calibration
CallBLKVOLA = blkimpv(Fwd, FiltOptTab.STRIKES(1:2:end), ZeroRate, TTM, FiltOptTab.MID(1:2:end), "Class", "call");
PutBLKVOLA  = blkimpv(Fwd, FiltOptTab.STRIKES(2:2:end), ZeroRate, TTM, FiltOptTab.MID(2:2:end), "Class", "put");

% plot of the smiles
% figure()
% plot(FiltOptTab.STRIKES(1:2:end), CallBLKVOLA, '*', 'DisplayName', 'Call Smile')
% grid on
% legend
% 
% figure()
% plot(FiltOptTab.STRIKES(2:2:end), PutBLKVOLA, '*', 'DisplayName', 'Put Smile')
% grid on
% legend

%% Call options

% delta between 0.1 and 0.9
d1             = 1./(CallBLKVOLA.*sqrt(TTM)).*log(Fwd./FiltOptTab.STRIKES(1:2:end)) + 0.5.*CallBLKVOLA.*sqrt(TTM);
CallDelta      = normcdf(d1); % verify if the discount is required
CallOneDateTab = FiltOptTab(1:2:end, :);
CallOneDateTab.BLKIMPVOLA = CallBLKVOLA;
OptIdx         = find((0.1 <= CallDelta).*(CallDelta <= 0.9));
CallOneDateTab = CallOneDateTab(OptIdx, :);

% OTM options
CallOneDateTab = CallOneDateTab(find(Fwd < CallOneDateTab.STRIKES), :);


%% Put Options

PutDelta      = abs(normcdf(d1) - 1);
PutOneDateTab = FiltOptTab(2:2:end, :);
PutOneDateTab.BLKIMPVOLA = PutBLKVOLA;
OptIdx        = find((0.1 <= PutDelta).*(PutDelta <= 0.9));
PutOneDateTab = PutOneDateTab(OptIdx, :);

% OTM options
PutOneDateTab = PutOneDateTab(find(Fwd > PutOneDateTab.STRIKES), :);

% output


end
