function y = nigpdf(x, delta, theta, k)
% 
% Function that gives the pdf of NIG distribution evaluated at input x
% 
% INPUT
% delta: average volatility param
% theta: skew param
% k:     vol of vol
% 
% OUTPUT
% y:     vector of pdf values at input points x
% 
% References:
% [1] Prause, K. (1999). The Generalized Hyperbolic Model
% 
% -------------------------------------------------------------------------

%% Constraints for the parameters
if delta <= 0
     error('delta must be positive.');
end
if k <= 0
     error('k must be positive.');
end


% transform input into column vector
[nx, mx] = size(x);
x = reshape(x, nx * mx, 1);

% % calculate leading factor
% q = (delta * alpha / pi) * exp(delta * sqrt(alpha^2 - beta^2));
% xbar = x - mu;
% z = sqrt(delta^2 + xbar .* xbar);
% 
% % modified Bessel function of the third kind 
% K = besselk(1, (alpha * z));
% 
% y = q ./ z .* K .* exp(beta * (x - mu));

% PDF calculation
BesselFun =  ...
    besselk(1, (sqrt(theta^2 + delta^2/k)/delta^2)...
        .*sqrt(x.^2 + delta^2/k));

y = (1/pi)*exp(1/k)*sqrt(theta^2/(k*delta^2) + 1/(k^2))*...
    exp((theta/(delta^2)).*x).*BesselFun./(sqrt((x.^2 + delta^2/k)));

% transform input vector back to input format
y = reshape(y, nx, mx);
y(isinf(x)) = 0;
y(isnan(y)) = 0;

end
