% Msc Thesis Masatti Gabriele
%
% Based on Patrizia Semeraro:
% Multivariate tempered stable additive subordination for financial models
%
% Supervisors: Michele Azzone, Roberto Baviera
% 
% -*- -*- -*- -*- -*- -*- -*- -*- -*-
%
% Academic Year 2022 - 2023
% Computation of the correlations between pair of assets
% Save the output in a .mat file
%

%% clearing
clear
close all
clc
rng(0);                                                                                 % We fix a seed to get my results replicable in case of randomness
warning('off');                                                                         % No warnings
options = optimset('Display','none');                                                   % Suppress function messages
set(0,'DefaultFigureWindowStyle','docked');                                             % Docking figures (more manageable)
format long                                                                             % More representation precision
FntNm = 'Times';
FntSz = 20;
formatData = 'dd/mm/yy'; 
IBDaycount = 3;

%% Correlations
AssetNames   = {'.SPX'; '.SX5E'};
TimeHorizons = [1; 5; 30; 60; 120];
Prices = xlsread("SPXN225Prices.xlsx", "PricesFilteredTrue");
Prices = Prices(:,2:end);
[Correlations] = AssetCorrs(Prices, TimeHorizons, AssetNames');

figure();
ax = gca;
plot(TimeHorizons, 100*Correlations(:, 2), '*-', 'DisplayName', 'Hist Correlation')
legend('FontName', FntNm, 'FontSize', FntSz)
grid on
title(strcat(AssetNames{1,1}, {' '}, AssetNames{2,1}, {' '}, 'Historical correlation vs TimeHorizon (in days)'), 'FontName', FntNm, 'FontSize', FntSz)
xlabel('TimeHorizon', 'FontName', FntNm, 'FontSize', FntSz)
ylabel('Correlation', 'FontName', FntNm, 'FontSize', FntSz)
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.YAxis.TickLabelFormat = '%g%%';
ax.XAxis.FontName = FntNm;
ax.YAxis.FontName = FntNm;

