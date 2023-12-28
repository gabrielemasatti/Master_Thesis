function [f1, f2] = simulateLogFwd(Nsim, M, flagspline, t, s, MarginalParams, CommonParameters)
%
% Function that simulates an increment of the log fwd process of the model 
% S-IG, using the common params calibrated on the historical correlations
%
% INPUT
% Nsim:             N of simulations required
% M:                approx of the FFT-S algorithm
% flagspline:       interpolation method in FFT
% t:                yearfrac up to the last time instant       
% s:                yearfrac up to the first time instant   
% MarginalParams:   [a_j, beta_j, delta_j] for both the ULs
% CommonParameters: [rho, a, q] common parameters
%
% OUTPUT
% f1: log-forward increment of asset 1
% f2: log-forward increment of asset 2
%

%% Preliminaries
A1a     = MarginalParams(1,1);
A1beta  = MarginalParams(1,2);
A1delta = MarginalParams(1,3);

A2a     = MarginalParams(2,1);
A2beta  = MarginalParams(2,2);
A2delta = MarginalParams(2,3);

rho = CommonParameters(1);
a   = CommonParameters(2);
q   = CommonParameters(3);

%% Simulation

% simulation of X1, X2 increment
X1 = simulateSATOIGFFT(Nsim, M, flagspline, t, s, A1a, A1beta, A1delta, a, q, 0);
X2 = simulateSATOIGFFT(Nsim, M, flagspline, t, s, A2a, A2beta, A2delta, a, q, 0);

% trying in a different way
X1lambda = 2*pi*(1-a*sqrt(A1a))^2;
X1mu     = sqrt(pi*A1a)*(1-a*sqrt(A1a));
X1distr  = (t^q-s^q)*random('inversegaussian', X1mu, X1lambda, Nsim, 1);

X2lambda = 2*pi*(1-a*sqrt(A2a))^2;
X2mu     = sqrt(pi*A2a)*(1-a*sqrt(A2a));
X2distr  = (t^q-s^q)*random('inversegaussian', X2mu, X2lambda, Nsim, 1);

% figure
% ksdensity(X1distr)
% hold on
% ksdensity(X1)
% grid on
% hold off
% 
% figure
% ksdensity(X2distr)
% hold on
% ksdensity(X2)
% grid on

% simulation of Z increment
Z       = simulateSATOIGFFT(Nsim, M, flagspline, t, s, A1a, A1beta, A1delta, a, q, 1);
Zlambda = 2*pi*a^2;
Zmu     = sqrt(pi)*a;
Zdistr  = (t^q-s^q)*random('inversegaussian', Zmu, Zlambda, Nsim, 1);

% figure
% ksdensity(Zdistr)
% hold on
% ksdensity(Z)
% grid on

% simulation of the indep Gaussians
G1 = randn(Nsim, 1);
G2 = randn(Nsim, 1);

% simulation of the corr Gaussians
Epsilon = [1, rho*sqrt(A1a)*sqrt(A2a)*A1delta*sqrt(A2delta); ...
    rho*sqrt(A1a)*sqrt(A2a)*A1delta*sqrt(A2delta), 1];

A = chol(Epsilon)';

G = (A*randn(2, Nsim))';

%% Computation of the increment
alpha = 0.5;  % generalize by giving as input
% c1t = gamma(-alpha).*((1/A1a)^alpha-(1/A1a - t.^q*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha);
% c2t = gamma(-alpha).*((1/A2a)^alpha-(1/A2a - t.^q*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha);
% c1s = gamma(-alpha).*((1/A1a)^alpha-(1/A1a - s.^q*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha);
% c2s = gamma(-alpha).*((1/A2a)^alpha-(1/A2a - s.^q*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha);

c1t = (gamma(-alpha).*((1/A1a)^alpha-(1/A1a - t.^q*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha).*(1-a.*sqrt(A1a)))+...
    (gamma(-alpha).*((1-(1 - t.^q*A1a*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha).*a));

c2t = (gamma(-alpha).*((1/A2a)^alpha-(1/A2a - t.^q*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha).*(1-a.*sqrt(A2a)))+...
    (gamma(-alpha).*((1-(1 - t.^q*A2a*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha).*a));


c1s = (gamma(-alpha).*((1/A1a)^alpha-(1/A1a - s.^q*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha).*(1-a.*sqrt(A1a)))+...
    (gamma(-alpha).*((1-(1 - s.^q*A1a*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha).*a));

c2s = (gamma(-alpha).*((1/A2a)^alpha-(1/A2a - s.^q*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha).*(1-a.*sqrt(A2a)))+...
    (gamma(-alpha).*((1-(1 - s.^q*A2a*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha).*a));


% f1 = A1beta*A1delta^2.*X1distr + A1delta.*sqrt(X1distr).*G1...
%     + A1beta*A1delta^2*A1a.*Zdistr + sqrt(A1a*A1delta^2).*sqrt(Zdistr).*G(:,1)+ (c1t-c1s);
% 
% f2 = A2beta*A2delta^2.*X2distr + A2delta.*sqrt(X2distr).*G2...
%     + A2beta*A2delta^2*A2a.*Zdistr + sqrt(A2a*A2delta^2).*sqrt(Zdistr).*G(:,2)+ (c2t-c2s);

f1 = A1beta*A1delta^2.*X1 + A1delta.*sqrt(X1).*G1...
    + A1beta*A1delta^2*A1a.*Z + sqrt(A1a*A1delta^2).*sqrt(Z).*G(:,1) + (c1t - c1s);

f2 = A2beta*A2delta^2.*X2 + A2delta.*sqrt(X2).*G2...
    + A2beta*A2delta^2*A2a.*Z + sqrt(A2a*A2delta^2).*sqrt(Z).*G(:,2) + (c2t - c2s);

%% check on the char function

A1phi_t  = @(z) exp(1i*z*c1t).*exp(2*sqrt(pi)*(sqrt(1/A1a) - sqrt(1/A1a - t^q*(1i*A1beta*A1delta^2*z - 0.5*A1delta^2*z.^2))));
A1phi_s  = @(z) exp(1i*z*c1s).*exp(2*sqrt(pi)*(sqrt(1/A1a) - sqrt(1/A1a - s^q*(1i*A1beta*A1delta^2*z - 0.5*A1delta^2*z.^2))));
A1phi_s_t = @(z) A1phi_t(z)./A1phi_s(z);
A1checkMart = A1phi_s_t(-1i);

A2phi_t  = @(z) exp(1i*z*c1t).*exp(2*sqrt(pi)*(sqrt(1/A2a) - sqrt(1/A2a - t^q*(1i*A2beta*A2delta^2*z - 0.5*A2delta^2*z.^2))));
A2phi_s  = @(z) exp(1i*z*c1s).*exp(2*sqrt(pi)*(sqrt(1/A2a) - sqrt(1/A2a - s^q*(1i*A2beta*A2delta^2*z - 0.5*A2delta^2*z.^2))));
A2phi_s_t = @(z) A2phi_t(z)./A2phi_s(z);
A2checkMart = A2phi_s_t(-1i);





end