function [X1, X2] = simulateSATONIG(Nsim, TTM, MarginalParams, CommonParameters, Method)
%
% Function that simulates the paths of S-IG Subordinated BM:
% It is possible to choose wether the paths of the underlyings should be
% correlated (params calibrated from hist correlations) of not.
%
% INPUT
% Nsim:             N of simulations required
% TTM:              Time up to which the simulations are required
% MarginalParams:   [a, beta, delta] in M(2x3) for the assets
% CommonParameters: [rho, a, q]
% Method:           1 if correlated, 0 if not (5.9) in Semeraro
%
% OUTPUT
% X1: log-increments of first asset
% X2: log-increments of second asset
%
%

switch Method

    case 0
       
        [X1, X2] = simulateDisjointSATONIG(Nsim, TTM, MarginalParams, CommonParameters(3));

    case 1

        [X1, X2] = simulateJointSATONIG(Nsim, TTM, MarginalParams, CommonParameters, Method);
        
end

