function [X1sato, X2sato] = simulateNIGCOPULA(Nsim, t, s, MarginalParams, rho, q)
%
% function that simulates the increments of two correlated assets with a
% Gaussian Copula, where the marginals are SSDNIG processes
%
% INPUT
% Nsim:           number of simulations required
% TTM:            yearfraction of the simulation
% MarginalParams: parameters of marginals [delta, beta, k]
% rho:            correlation parameter of the gaussian copula
%
% OUTPUT
% X1sato: increment of the logreturn of first asset
% X2sato: increment of the logreturn of second asset
%
%


%% Preliminaries
TTM     = t-s;
A1delta = MarginalParams(1,1);
A1theta = MarginalParams(1,2);
A1k     = MarginalParams(1,3);

A2delta = MarginalParams(2,1);
A2theta = MarginalParams(2,2);
A2k     = MarginalParams(2,3);

%% Simulate a 2-dim gaussian vector with correlation rho
rhoMatrix = eye(2,2) + rho*(triu(ones(2,2)) - eye(2,2)) + rho*(tril(ones(2,2)) - eye(2,2));
A  = chol(rhoMatrix)';

% correlated gaussians
G = (A*randn(2, Nsim))';

% map to u1 and u2
U = normcdf(G);

% risk-free corrections
charFun1_t = @(u) exp(1/A1k *(1 - sqrt(1 + u.^2.*t^(2*q)*A1delta^2*A1k - 2i.*A1theta.*t^q.*u*A1k)));
charFun1_s = @(u) exp(1/A1k *(1 - sqrt(1 + u.^2.*s^(2*q)*A1delta^2*A1k - 2i.*A1theta.*s^q.*u*A1k)));
charFun1_s_t = @(u) charFun1_t(u)./charFun1_s(u);


charFun2_t = @(u) exp(1/A2k *(1 - sqrt(1 + u.^2.*t^(2*q)*A2delta^2*A2k - 2i.*A2theta.*t^q.*u*A2k)));
charFun2_s = @(u) exp(1/A2k *(1 - sqrt(1 + u.^2.*s^(2*q)*A2delta^2*A2k - 2i.*A2theta.*s^q.*u*A2k)));
charFun2_s_t = @(u) charFun2_t(u)./charFun2_s(u);

c1_t = -log(charFun1_t(-1i));
c1_s = -log(charFun1_s(-1i));
c1_s_t = -log(charFun1_s_t(-1i));

c2_t = -log(charFun2_t(-1i));
c2_s = -log(charFun2_s(-1i));
c2_s_t = -log(charFun2_s_t(-1i));

% invert the NIG cdf
X1sato = simulateATSFFT(15, Nsim, 1, t, s, A1delta, A1delta, A1theta, A1theta, ...
    A1k, A1k, U(:, 1), q) + c1_s_t;

X2sato = simulateATSFFT(15, Nsim, 1, t, s, A2delta, A2delta, A2theta, A2theta, ...
    A2k, A2k, U(:, 2), q) + c2_s_t;

end