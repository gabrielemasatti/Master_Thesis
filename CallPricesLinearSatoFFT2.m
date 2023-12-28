function [prices] = CallPricesLinearSatoFFT2(forward, discount, moneyness, timeToMaturity, param, numericalMethodParameters, alpha, nAsset)
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

        C1 = param(1);
        C2 = param(2);
        lambda1 = param(3);
        lambda2 = param(4);

    case 2

        C1 = param(5);
        C2 = param(6);
        lambda1 = param(7);
        lambda2 = param(8);

end

q = param(9);

% char function of X_j(t)
charFunY = @(u) exp(-2*sqrt(pi).*(C1.*(sqrt(lambda1 - 1i.*timeToMaturity^q.*u) - sqrt(lambda1))))...
    .* exp(-2*sqrt(pi).*(C2.*(sqrt(lambda2 + 1i.*timeToMaturity^q.*u) - sqrt(lambda2))));

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


