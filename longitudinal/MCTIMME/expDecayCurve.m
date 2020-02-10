function [v] = expDecayCurve(p,times)
% calculate exponential decay curve, used to estimate relaxation time
% hyperparameters
%
% inputs:
% p = vector of function parameters
% times = vector of time points to evaluate function at
%
% output:
% v = vector of values at which function has been evaluated
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% expDecayCurve.m (version 1.00)

v = p(2)*exp(-p(1)*times) + (1-exp(-p(1)*times))*p(3);
v = v';

