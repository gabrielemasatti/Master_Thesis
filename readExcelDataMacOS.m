function [optionData] = readExcelDataMacOS(filename, sheetname)
% 
% Reads data from excel (MacOS version)
% It reads data from excel
%
% INPUT
% filename: excel file name where data are stored
% 
% OUTPUTS:
% OptionStruct: struct with settlement, expiries, ask price, bid price,
%               strike price, open interest, type of option
%

%% Dates from Excel

% Define cell-to-array function
cell2arr = @(c) reshape(c, size(c));

% Import all content
allcontent = readcell(filename,"Sheet", sheetname);

% Settlement date
settlement            = allcontent{2,7};            
optionData.settlement = datenum(settlement);

% Maturity dates of the options
maturities   = allcontent(3:end, 6);
nOptions     = size(maturities, 1);
maturityVec  = zeros(nOptions, 1);
for i = 1:nOptions
    maturityVec(i) = datenum(matlab.datetime.compatibility.convertDatenum(maturities{i,1}));
end
optionData.maturities = maturityVec;

% Bid-Ask-Opening prices
bidPrices            = cell2arr(allcontent(3:end, 3));
askPrices            = cell2arr(allcontent(3:end, 4));
strikes              = cell2arr(allcontent(3:end, 5));
bidPrices            = str2double(bidPrices);
askPrices            = str2double(askPrices);
strikes              = str2double(strikes);

idxBidPrices = isnan(bidPrices);
idxaskPrices = isnan(askPrices);
idxstrikes   = isnan(strikes);

bidPrices(idxBidPrices) = 0;
askPrices(idxaskPrices) = 0;
strikes(idxstrikes)     = 0;

optionData.bidPrices = bidPrices;
optionData.askPrices = askPrices;
optionData.strikes   = strikes;


% Call or put
% flag                 = cell2arr(allcontent(3:end, 13));
% optionData.flag      = flag;

end % readExcelDataOSX