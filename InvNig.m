function [ x,error] = InvNig(y,pars,dx,n,tol )
% finds the inverse of NIG distribution using Newton method
% y --->> the input the inverse cdf of y is the objective
% x0 --->> initial guess
% pars = [alpha,beta,mu,delata,Xmin,Xmax]
% n -->> maximum number of iteration
% tol -->> the tolerance for the error part 

% email:        iasadzad@ucalgary.ca    
% 
% Implementation Date:  15 May,2015
% Author:              Ilnaz Asadzadeh
% -------------------------------------------------
if nargin ==3
    tol = 1e-4;
    n   = 1e1;
end
alpha    = pars(1); beta = pars(2); mu = pars(3); delta = pars(4);
Xmin     = pars(5); Xmax = pars(6);
Phi0     = NIGcdf(0,alpha,beta,mu,delta);
f        = @(x) NIGcdf(x,alpha,beta,mu,delta,Phi0)- y;
g        = @(x) x - f(x)/nigpdf(x,alpha,beta,mu,delta);
% first try
tryagain = 1;
while tryagain
    [x,bracket]   = invNIGinitial(y,Xmin,Xmax,tol,pars,dx);
    testbracket   = @(x) ( (x-bracket(1))*(bracket(2)-x)<0 );
    tryagain      = 0;
    error(1)      = min(abs(x-bracket));
    k = 1;
    while ~tryagain && error(k) >= tol && k<= n
        %disp(k)
        x(k+1)      = g(x(k));
        error(k+1) = abs(x(k+1)-x(k));
        if testbracket(x(k+1))
            tryagain = 1;
            dx = dx/2; 
            Xmin = bracket(1);
            Xmax = bracket(2);
            %disp(dx);
            %keyboard
            % the loop will end and we will re-initialize
        end
        k          = k+1;
    end
end
x = x(end); error = error(end);

end

