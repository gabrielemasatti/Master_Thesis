close all
clear 
warning off

% parameters 
A1delta = 0.1329;
A1theta = -0.1034;
A1k     = 2.0831;

A2delta = 0.1485;
A2theta = -0.1021;
A2k     = 1.8819;

q = 0.688;
TTM = 3/12;
% rho = 0.855; 1Y param
rho = 0.81;
xx = linspace(-2, 2, 50);

%% Build z1 and z2 integration extrema as in (3.4) of paper

A1pdf = @(x) nigpdf(x, A1delta, A1theta, A1k);
A1pdf = @(x) A1pdf(x).*(A1pdf(x) >= 1e-4);
% A1pdf = A1pdf(xx);
% figure
% plot(xx, A1pdf)

A2pdf = @(x) nigpdf(x, A2delta, A2theta, A2k);
A2pdf = @(x) A2pdf(x).*(A2pdf(x) >= 1e-4);
% A2pdf = A2pdf(xx);
% hold on
% plot(xx, A2pdf)
% hold off

A1cdf      = @(x1) nigcdf(x1, A1delta, A1theta, A1k);
A1cdf    = @(x1) A1cdf(x1).*(A1cdf(x1) >= 0).*(A1cdf(x1) <= 1) +...
    0.999999999999.*(A1cdf(x1) >= 1);
A1cdfres = A1cdf(xx);
figure
plot(xx, A1cdfres, 'DisplayName', 'x1')
hold on

A2cdf      = @(x2) nigcdf(x2, A2delta, A2theta, A2k);
A2cdf    = @(x2) A2cdf(x2).*(A2cdf(x2) >= 0).*(A2cdf(x2) <= 1) +...
    0.999999999999.*(A2cdf(x2) >= 1);
A2cdfres = A2cdf(xx);
plot(xx, A2cdfres, 'DisplayName', 'x2')
hold on

%% Build the joint cdf
A1z    = @(x1) norminv(A1cdf(x1));
A1zres = A1z(xx);

figure
plot(xx, A1zres, 'DisplayName', 'x1')
hold on

A2z  = @(x2) norminv(A2cdf(x2));
A2zres = A2z(xx);

plot(xx, A2zres, 'DisplayName', 'x2')
grid on
legend
hold off

GaussPDF1 = @(x1) 1./(normpdf(A1z(x1)));
GaussPDF1res = GaussPDF1(xx);
GaussPDF2 = @(x2) 1./(normpdf(A2z(x2)));
GaussPDF2res = GaussPDF2(xx);

MultiPDF  = @(x1, x2) 1./(2*pi*sqrt(1 - rho^2)).*...
    exp(-(A1z(x1).^2 - 2*rho*A1z(x1).*A2z(x2) + A2z(x2).^2)./(2*(1-rho^2)));

NIGPdf1 = @(x1) A1pdf(x1);
NIGPdf2 = @(x2) A2pdf(x2);

JointPDF = @(x1, x2) MultiPDF(x1, x2)...
    .*NIGPdf1(x1).*NIGPdf2(x2).*GaussPDF1(x1).*GaussPDF2(x2);

[X, Y] = meshgrid(xx, xx);
figure
surf(xx, xx, JointPDF(X, Y))

%%% TRUNCATION DOMAIN %%%
% Use as truncation for the integral CDF 

%% Correlation Function
% CDFcheck  = integral2(@(x1,x2) JointPDF(x1, x2), -2.9, 2.6, -2.9, 2.6) % 1Y
CDFcheck  = integral2(@(x1,x2) JointPDF(x1, x2), -2.8, 2.8, -2.8, 2.8);

% CDFcheck  = @(x1) quadgk(@(x2) JointPDF(x1, x2), -2, 2);
% CDFcheckB = @(x1) arrayfun(CDFcheck, x1);
% CDFcheck  = quadgk(CDFcheckB, -20, 2);

% covariance
% EX1X2_1 = integral2(@(x1,x2) x1.*x2.*JointPDF(x1, x2), -2.8, 2.5, -2.8, 2.5); %1Y
EX1X2_1 = integral2(@(x1,x2) x1.*x2.*JointPDF(x1, x2), -2.8, 2.8, -2.8, 2.8); %1Y

% EX1X2_1_check = TTM^(4*q).*integral2(@(x1,x2) x1.*x2.*JointPDF(x1.*TTM^q, x2.*TTM^q), -2.5, 2.5, -2.5, 2.5); 
% EX1X2_1 = @(x1) quadgk(@(x2) x1.*x2.*JointPDF(x1, x2), -0.2, 0.2);
% B       = @(x1) arrayfun(EX1X2_1, x1);
% EX1X2_2 = quadgk(B, -20, 0.7);

% other moments
EX1 = A1theta;
EX2 = A2theta;
V1  = A1delta^2 + A1theta^2*A1k;
V2  = A2delta^2 + A2theta^2*A2k;

% Corr Function
rhoFun = (EX1X2_1 - EX1*EX1)/sqrt(V1*V2);

