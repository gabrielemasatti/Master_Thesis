function [ModelCorrelations] = ComputeCorrelation(Asset1Params, Asset2Params, FracHoriz, rho, a, q)
%
% Function that computes the model correlations on different time horizons
% 
% INPUT
% Asset1Params: Marginal parameters of asset1 [a1, beta1, delta1]
% Asset2Params: Marginal parameters of asset2 [a2, beta2, delta2]
% FracHoriz:    Time Horizons on which the model corr is computed
% rho:          rho_{ij} common parameter
% a:            common parameter of the common Z subordinator
% q:            common self-similar parameter
%
% OUTPUT
% ModelCorrelations: model correlations on different time horizons
% 

% moments of the RVs considered
Zmean  = sqrt(pi)*a;
ZVar   = 0.5*sqrt(pi)*a;
S1mean = sqrt(pi*Asset1Params(1));
S1Var  = 0.5*sqrt(pi*Asset1Params(1)^3);
S2mean = sqrt(pi*Asset2Params(1));
S2Var  = 0.5*sqrt(pi*Asset2Params(1)^3);

ModelCorrelations = zeros(size(FracHoriz,1), 1);

% computation of the correlation according to (1.16) of paper
for i = 1:size(ModelCorrelations,1)

    ModelCorrelations(i) = (rho.*Asset2Params(3).*Asset1Params(3).*sqrt(Asset2Params(1).*Asset1Params(1)).*Zmean.*FracHoriz(i).^q+ ...
        Asset2Params(2).*Asset1Params(2).*Asset2Params(3).^2.*Asset1Params(3).^2.*Asset2Params(1).*Asset1Params(1).*FracHoriz(i).^(2*q).*ZVar)...
        ./(sqrt((Asset2Params(3).^2.*FracHoriz(i).^q*S2mean+ Asset2Params(3).^4*Asset2Params(2)^2*S2Var.*FracHoriz(i).^(2*q))...
        *(Asset1Params(3)^2.*FracHoriz(i).^q*S1mean+ Asset1Params(3)^4*Asset1Params(2)^2*S1Var.*FracHoriz(i).^(2*q))));

end

end