function [yCDF1] = NIGCDFComputation(MarginalParams, q, x)
%
% Function that returns the CDF of a NIG
%
% INPUT
% MarginalParams: MarginalParams(1) -> delta
%                 MarginalParams(2) -> theta
%                 MarginalParams(3) -> k
% q: self-similar parameter
% 
% OUTPUT
% yCDF: cdf evaluated in the x point
%

% NIG Parameters
delta = MarginalParams(1);
theta = MarginalParams(2);
k     = MarginalParams(3);

% NIG pdf according to Cont-Tankov
BesselFun = @(x) ...
    besselk(1, (sqrt(theta^2 + delta^2/k)/delta^2)...
        .*sqrt(x.^2 + delta^2/k));

pdf = @(x) (1/pi)*exp(1/k)*sqrt(theta^2/(k*delta^2) + 1/(k^2))*...
    exp((theta/(delta^2)).*x).*BesselFun(x)./(sqrt((x.^2+delta^2/k)));

% NIG CDF
res   = @(y) quadgk(pdf, -10, y);
yCDF1 = res(x);

end