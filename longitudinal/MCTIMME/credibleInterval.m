classdef credibleInterval
% object to define a 95% credible interval
%
% Microbial Counts Trajectories Infinite Mixture Model Engine (MC-TIMME)
% Copyright (C) 2012 Georg Gerber
% credibleInterval.m (version 1.00)
%
    properties
        med % median
        q025 % 2.5% quantile
        q975 % 97.5% quantile
    end
end

