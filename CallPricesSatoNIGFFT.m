function [prices] = CallPricesSatoNIGFFT( forward, discount, moneyness, timeToMaturity, params, numericalMethodParameters, alpha, q)
%
% Compute the prices of Calls with FFT method
%
% INPUT 
%
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
delta = params(1);
theta = params(2);
k     = params(3);
% q     = params(4);

% char function of Y^rho_j
% charFunY = @(u) exp(mu*(1 - sqrt(1 - 2i*(u.*theta + 0.5*delta^2*u.^2*1i).*timeToMaturity^q)));
charFunY = @(u) exp(1/k *(1 - sqrt(1 + u.^2.*timeToMaturity^(2*q)*delta^2*k- 2i.*theta.*timeToMaturity^q.*u*k)));

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