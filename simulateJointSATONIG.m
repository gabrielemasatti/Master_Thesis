function [X1, X2] = simulateJointSATONIG(Nsim, t, MarginalParams, CommonParams)
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

rho = CommonParams(1);
a   = CommonParams(2);
q   = CommonParams(3);

%% Simulation of the IG RVS with inverse cdf




end