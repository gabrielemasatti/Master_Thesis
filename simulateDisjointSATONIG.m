function [X1, X2] = simulateDisjointSATONIG(Nsim, t, MarginalParams, q)
%
% Function that simulates an S - IG increment Nsim times
% 
% INPUT
% Nsim:           N of simulations required
% t:              Time up to which the simulations are required
% MarginalParams: [a, beta, delta] in M(2x3) for the assets
% q:              self-similar parameter          
%
% OUTPUT
% X1: log-increments of first asset
% X2: log-increments of second asset
%%

%% Preliminaries
A1a     = MarginalParams(1,1);
A1beta  = MarginalParams(1,2);
A1delta = MarginalParams(1,3);

A2a     = MarginalParams(2,1);
A2beta  = MarginalParams(2,2);
A2delta = MarginalParams(2,3);

%% Simulation of the IG
% NIG params relationships
% A1lambda = 2*pi*(1-a*sqrt(A1a))^2;
% A2lambda = 2*pi*(1-a*sqrt(A2a))^2;
% A1mu = 1;
% A2mu = 1;
% A1mu = sqrt(A1a)*sqrt(pi)*(1-a*sqrt(A1a));
% A2mu = sqrt(A2a)*sqrt(pi)*(1-a*sqrt(A2a));
% A1lambda = 2*pi*(1-a*sqrt(A1a))^2;
% A2lambda = sqrt(pi)*(1-a*sqrt(A2a))*sqrt(A2a);
% A1mu = 1;
% A2mu = 1;

% simulate indep RVs needed 
% U1 = unifrnd(0, 1, [Nsim, 1]); % UNIF(0, 1)
% U2 = unifrnd(0, 1, [Nsim, 1]); % UNIF(0, 1)
% N1 = randn(Nsim, 1);
% N2 = randn(Nsim, 1);
% Y1 = N1.^2;
% Y2 = N2.^2;
% 
% L1 = A1mu + 0.5 .* A1mu.^2 .* Y1 ./ A1lambda - 0.5 * A1mu ./ A1lambda * sqrt( 4 .* A1mu .* A1lambda .* Y1 + A1mu.^2 * Y1 .^2 );
% A1temp = A1mu ./ ( L1 + A1mu );
% 
% L2 = A2mu + 0.5 .* A2mu.^2 .* Y2 ./ A2lambda - 0.5 * A2mu ./ A2lambda * sqrt( 4 .* A2mu .* A2lambda .* Y2 + A2mu.^2 * Y2 .^2 );
% A2temp = A2mu ./ ( L2 + A2mu );
% 
% A1indexes1 = U1 <= A1temp;
% A1indexes2 = U1 > A1temp;

% X1j             = random('inversegaussian', A1lambda, A1mu, Nsim, 1);
% X1j(A1indexes1) = L1(A1indexes1);
% X1j(A1indexes2) = A1mu.^2 ./ L1(A1indexes2);
% X1j = t^q*X1j;
% 
% A2indexes1 = U2 <= A2temp;
% A2indexes2 = U2 > A2temp;

% X2j             = random('inversegaussian', A2lambda, A2mu, Nsim, 1);
% X2j(A2indexes1) = L2(A2indexes1);
% X2j(A2indexes2) = A1mu.^2 ./ L2(A2indexes2);
% X2j = t^q*X2j;

% Z
% Zlambda = 2*sqrt(pi);
% Zmu = 1;
% Zlambda = 2*pi*a^2;
% Zmu = sqrt(pi)*a;
% 
% ZU = unifrnd(0, 1, [Nsim, 1]); % UNIF(0, 1)
% ZN = randn(Nsim, 1);
% ZY = ZN.^2;
% 
% ZL = Zmu + 0.5 .* Zmu.^2 .* ZY ./ Zlambda - 0.5 * Zmu ./ Zlambda * sqrt( 4 .* Zmu .* Zlambda .* ZY + Zmu.^2 * ZY .^2 );
% Ztemp = Zmu ./ ( ZL + Zmu );
% 
% Zindexes1 = ZU <= Ztemp;
% Zindexes2 = ZU > Ztemp;

% Z            = random('inversegaussian', Zlambda, Zmu, Nsim, 1);
% Z(Zindexes1) = ZL(Zindexes1);
% Z(Zindexes2) = Zmu.^2 ./ ZL(Zindexes2);
% Z = t^q.*Z;

S1mu     = 2*pi;
S1lambda = sqrt(pi*A1a);
S1       = random('inversegaussian', S1lambda, S1mu, Nsim, 1);
S1 = t^q*S1;

S2mu     = 2*pi;
S2lambda = sqrt(pi*A2a);
S2       = random('inversegaussian', S2lambda, S2mu, Nsim, 1);
S2 = t^q*S2;

% figure
% ksdensity(X1j)
% hold on 
% ksdensity(X2j)
% hold on
% ksdensity(Z)
% hold off

%% Simulation of the Gaussians and computation of the increment
% rhoMatrix = eye(2,2) + rho*(triu(ones(2,2)) - eye(2,2)) + rho*(tril(ones(2,2)) - eye(2,2));
% Epsilon = [1 rho*A1delta*A2delta*sqrt(A1a)*sqrt(A2a);...
%     rho*A1delta*A2delta*sqrt(A1a)*sqrt(A2a), 1];
% A  = chol(Epsilon)';

% independent gaussians
G1 = randn(Nsim, 1);
G2 = randn(Nsim, 1);

% correlated gaussians
% G = (A*randn(2, Nsim))';

% correct drift
alpha = 0.5;  % generalize by giving as input
c1 = gamma(-alpha).*((1/A1a)^alpha-(1/A1a - t.^q*(A1beta*A1delta^2 + 0.5*A1delta^2)).^alpha);
c2 = gamma(-alpha).*((1/A2a)^alpha-(1/A2a - t.^q*(A2beta*A2delta^2 + 0.5*A2delta^2)).^alpha);

%% Output
% X1 = (A1beta*A1delta^2*X1j + A1delta*sqrt(X1j).*G1)...     % uncorr part
%     + (A1beta*A1delta^2*A1a*Z + sqrt(A1a*A1delta^2)*sqrt(Z).*G(:,1)) + c1; % corr part
% 
% X2 = (A2beta*A2delta^2*X2j + A2delta*sqrt(X2j).*G2)...     % uncorr part
%     + (A2beta*A2delta^2*A2a*Z + sqrt(A2a*A2delta^2)*sqrt(Z).*G(:,2)) + c2; % corr part


X1 = (A1beta*A1delta^2*S1 + A1delta*sqrt(S1).*G1)+c1;

X2 = (A2beta*A2delta^2*S2 + A2delta*sqrt(S2).*G2)+c2;

end