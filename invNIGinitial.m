function [x0,bracket] = invNIGinitial(y,Xmin,Xmax,tol,pars,dx,Phi0)
% Finds close value of the inverse of cumalative distribution function
% Nigcdf(x0) ~= y, y \in (0,1)
% Xmin -->> minimum value of X (data)
% Xmax -->> maximum value of X (data)
% pars -->> (alpha, beta, mu, delta)
% dx -->> step size for X grid 
if nargin < 4 
    tol = 10^(-4);
end
if y ==1;
    x0 = Xmax;
elseif y==0
    x0 = Xmin;
end 
alpha   = pars(1); beta = pars(2);  mu = pars(3) ; delta = pars(4);
X       = Xmin:dx:Xmax;
X = unique([X,Xmax]);
if nargin>6
Phi     = NIGcdf(X,alpha,beta,mu,delta,Phi0);
else
Phi     = NIGcdf(X,alpha,beta,mu,delta);
end
f       = Phi - y; j = find(f>0); j= j(1)-1 ; 
if f(j) >= -tol && f(j+1) <= tol
dispay Done
x0 = X(j);
else 
    lambda = f(j+1)/(X(j+1)-X(j));
    x0     = lambda*X(j)+ (1-lambda)*X(j+1);
end
bracket = [X(j),X(j+1)];
if (x0-X(j))*(X(j+1)-x0)<0, keyboard, end

end

