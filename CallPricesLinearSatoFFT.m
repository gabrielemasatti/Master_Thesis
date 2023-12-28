function [prices] = CallPricesLinearSatoFFT(forward, discount, moneyness, timeToMaturity, param, numericalMethodParameters, alpha, nAsset)
%
% Compute the prices of Calls with FFT method
%
% INPUT 
% forward:                   forward value at settelment date
% discount:                  discount at maturity
% moneyness:                 grid of the moneyness
% timeToMaturity:            time to maturity
% [sigma, eta, k]:           set of marginal parameters
% numericalMethodParameters: parameters of the numerical method as a struct x1, xN, dx, z1, zN, dz, M for fft
% nAsset:                    evaluate correct params according to case 1/2
% 
% OUTPUT
% prices: Call option prices
%
% CALLS
% computeIntegral
%

switch nAsset
    
    case 1

        sigma = param(1);
        eta = param(2);
        k = param(3);

    case 2

        sigma = param(4);
        eta = param(5);
        k = param(6);

end

q = param(7);

% char function of X_j(t)
charFunY = @(u) exp(1/k *(1 - sqrt(1 + u.^2.*timeToMaturity^(2*q)*sigma^2*k- 2i.*eta.*timeToMaturity^q.*u*k)));

% MARTINGALE CONDITION exponent c 
c = -log(charFunY(-1i));

% char function
phi = @(u) exp(1i.*u.*c).*charFunY(u);

% integrand function
fTS            = @(v) 1./(v.^2 + 0.25);
LevyFunctionTS = @(v) 1/(2*pi).*phi(-v-1i/2).*fTS(v);
Integrals      = computeIntegral(LevyFunctionTS, moneyness, [], numericalMethodParameters, 1);
prices         = ((1 - exp(-moneyness./2).* Integrals)*discount*forward)';

end


