function [prices] = CallPricesSatoFFT( forward, discount, moneyness, timeToMaturity, params, numericalMethodParameters, alpha, q)
%
% Compute the prices of Calls with FFT method
%
% INPUT 
% forward:                   forward value at settelment date
% discount:                  discount at maturity
% moneyness:                 grid of the moneyness
% timeToMaturity:            time to maturity
% params:                    a, beta, delta of SATO NIG process
% numericalMethodParameters: parameters of the numerical method as a struct x1, xN, dx, z1, zN, dz, M for fft
%
% OUTPUT
% prices: Call option prices
%
% CALLS
% computeIntegral
%

%% model params
a     = params(1);
beta  = params(2);
delta = params(3);

% char function of Y^rho_j
charFunY = @(u) exp(gamma(-alpha).*((1./a - timeToMaturity.^q.*(1i.*beta*delta^2.*u - 0.5.*delta.^2.*u.^2)).^alpha...
    - (1./a).^alpha));

% MARTINGALE CONDITION exponent c 
c = gamma(-alpha).*((1/a)^alpha-(1/a - timeToMaturity.^q*(beta*delta^2 + 0.5*delta^2)).^alpha);

% char function
phi = @(u) exp(1i.*u.*c).*charFunY(u);

% integrand function
fTS            = @(v) 1./(v.^2 + 0.25);
LevyFunctionTS = @(v) 1/(2*pi).*phi(-v-1i/2).*fTS(v);
Integrals      = computeIntegral(LevyFunctionTS, moneyness, [], numericalMethodParameters, 1);
prices         = ((1 - exp(-moneyness./2).* Integrals)*discount*forward)';

end


