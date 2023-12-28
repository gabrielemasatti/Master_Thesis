function [ P ] = NIGcdf(x,alpha,beta,mu,delta,Phi0 )
x         = x(:)';
N         = size(x,2);
P         = zeros(size(x));
pdfHandle = @(x) nigpdf(x,alpha,beta,mu,delta);

%P         = arrayfun(@(x) quadgk(pdfHandle,-Inf,x),x);
%if length(x)>1, keyboard,end
if nargin>5 % we have a pre-calculated value we can use
    P = [Phi0,P];
    x = [0,x];
    N = N+1;
else
    % we want to avoid an infinite domain if we can
    P(1) = quadgk(pdfHandle,-inf,x(1));
end
for k = 2 : N
    P(k) = P(k-1) + quadgk(pdfHandle,x(k-1),x(k));
end
if nargin>5
    P = P(2:end);
end


end

