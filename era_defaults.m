function era_prefs = era_defaults
%
%Default settings for the ERA Toolbox
%
%Last Updated 7/21/16
%

%This script will be read by the ERA Toolbox to define the default settings
%for analyzing and viewing data in the ERA Toolbox. Should the user want to
%change any default settings, the values below can be changed. Upon
%restarting the toolbox, these settings will be read. Comments are provided 
%for understanding each input.
%
%Note: These changes will be made to each era_defaults.m every time a new 
%version is downloaded. 
%
%Input
% No inputs required in the command line
%
%Output
% era_prefs - data structure containing 
%

% Copyright (C) 2016 Peter E. Clayson
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see 
%     <http://www.gnu.org/licenses/>.
%

%History
% by Peter Clayson (7/21/16)
% peter.clayson@gmail.com
%

%define an empty structure array where preferences will be stored
era_prefs = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Processing Preferences%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of chains for processing data in CmdStan (default = 3)
era_prefs.proc.nchains = 3;

%Number of iterations for processing data in CmdStan (default = 1000)
era_prefs.proc.niter = 1000;

%Verbose input while processing data in CmdStan (defaults = 1);
% 1 = No
% 2 = Yes
era_prefs.proc.verbose = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Viewing Preferences%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Value to be used for dependability cutoff/threshold
era_prefs.view.depvalue = .80;

%For the inputs below a value of 1 (default) indicates that the
%figure/table should be viewed; a value of 0 indicates the figure/table
%should not be viewed
%The default value for all of these figures/tables is 1

%Plot that shows the relationship between the number of trials retained 
%for averaging and dependability
era_prefs.view.plotdep = 1;

%Plot that compares the intraclass correlation coefficients for each
%group and/or condition
era_prefs.view.ploticc = 1;

%Table displaying information about cutoffs based on dependability
%threshold
era_prefs.view.inctrltable = 1;

%Table displaying information about overall dependability with data
%including all trials after applying the trial cutoffs
era_prefs.view.overalltable = 1;

%Table displaying information about between- and within-person standard
%deviations as well as intraclass correlation coefficients
era_prefs.view.showstddevt = 1;

%Figure displaying between-person standard deviations for each group and/or
%condition
era_prefs.view.showstddevf = 1;


%Which dependability estimate to plot on the figure that shows the
%relationship between the number of trials retained for averaging and
%dependability (default = 2)
% 1 = lower limit of the credible interval
% 2 = point estimate of the credible interval
% 3 = upper limit of the credible interval
era_prefs.view.plotdepline = 2;

%Number of trials to plot on the figure that shows the relationship between
%the number of trials retained for averaging and dependability 
%(default = 50)
era_prefs.view.ntrials = 50;

%Which cutoff to use for estimating the dependability of averages for a
%trial cutoff (default = 2);
% 1 = lower limit of the credible interval
% 2 = point estimate of the credible interval
% 3 = upper limit of the credible interval
era_prefs.view.meascutoff = 2;

%Measure of central tendency to use for estimating the overall score
%dependability of waveforms after applying the trial cutoffs (default = 1)
% 1 = mean
% 2 = median
era_prefs.view.depcentmeas = 1;

end