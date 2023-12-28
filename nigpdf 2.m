function y = nigpdf(x, alpha, beta, mu, delta)
%  NIGpdf Probability density function (PDF) for Normal-Inverse-Gaussian distribution.
%   y = NIGpdf(x, alpha, beta, mu, delta) returns the pdf of the 
%   Normal-Inverse-Gaussian distribution with the parameter beta which 
%   determines the skewness, alpha the shape parameter, mu the location parameter
%   and delta scale parameter, evaluated at the values in x.
%   
%   alpha, beta, mu, delta could be scalar or vector for time varying distribution.
%   alpha and delta must be positive values.
%   alpha > |beta| must hold.


%% The pdf
% transform input into column vector
x    = x(:)';
p    = (delta * alpha / pi) .* exp(delta * sqrt(alpha.^2 - beta.^2));
xbar = x - mu;
z    = sqrt(delta.^2 + xbar .* xbar);
K    = besselk(1, (alpha .* z));
y    = p ./ z .* K .* exp(beta .* (x - mu));

