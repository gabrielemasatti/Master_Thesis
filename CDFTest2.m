clear all
clc

q = 0.6; %(to be passed as input in common params)

% Marginal Parameters
A1Delta = 0.25;
A1Theta = -12;
A1k     = 10;

A2Delta = 0.10;
A2Theta = -30;
A2k     = 2;

xx = [0.001:0.01:3]';
yy = [0.001:0.01:3]';

% Build CDFs
FracHoriz  = 1./365;
PDF = @(y) 1./(sqrt(2*pi*A1k*y.^3)).*exp(-(y - 1).^2./(2*A1k.*y)).*(y>0);
IntFun = @(x, y) ...
    normcdf((x'./FracHoriz^q - A1Theta*y')./(A1Delta.*sqrt(y'))).*PDF(y)';
CDF = @(x) quadgk(@(y) IntFun(x, y)', 0.1, 1000);

yres = zeros(length(xx), 1);

for i=1:length(xx)
    yres(i) = CDF(xx(i));
end