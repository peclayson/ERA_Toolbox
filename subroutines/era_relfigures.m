function era_relfigures(varargin)
%Creates various figures and tables for dependability data. See user manual
% for more specific information about each figure.
%
%era_relfigures('era_data',era_data)
%
%Last Modified 8/17/17
%
%Required Inputs:
% era_data - ERA Toolbox data structure array containing outputs from
%  CmdStan
%
%Optional Inputs:
% depcutoff - dependability level to use for cutoff when deciding the
%  minimum number of trials needed to achieve this specified level of
%  dependability
% plotdep - plot the table displaying dependability v number of trials
%  included in average
% ploticc - plot the intraclass correlation coefficients for data. This ICC
%  can be interpreted as the proportion of total variance that is between
%  persons
% showinct - display table with information for each event/group
%  combination (cutoffs and dependability)
% showoverallt - display table with information for data including all
%  trials for those participants that meet cutoff threshhold
% showstddevt - display table with information about the sources of
%  variance (between person v within person)
% plotbetstddev - plot the between-person standard deviations stratified
%  by group and event
% plotdepline - indicate whether to plot 1-lower limit of credible
%  interval, 2-point estimate, 3-upper limit of credible interval for the  
%  plot of dependability v number of trials (default: 1)
% plotntrials - indicate the number of trials to plot (x-axis) in the 
%  dependability v number of trials plot (default: 50)
% meascutoff - which estimate to use to define cutoff for number of trials
%  1 - lower limit of credible interval, 2 - point estimate, 3 - upper
%  limit of credible interval (default: 1)
% depcentmeas - which measure of central tendency to use to estimate the
%  overall dependability, 1 - mean, 2 - median (default: 1)
%
%Output:
% Various figures may be plotted. Tables will also have buttons for
%  saving the table to a file. Additionally the table showing the cutoffs
%  for numbers of trials to achieve a given level of dependability will
%  have a button to save the ids of participants to include and exclude.

% Copyright (C) 2016-2018 Peter E. Clayson
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
% by Peter Clayson (4/18/16)
% peter.clayson@gmail.com
%
%4/20/16 PC
% bug fix: point estimate and lower limit were flipped for dependability
%  plot
% changed spacing in guis a bit
% add median number of trials to the overall inclusion table (gui and
%  output table)
% aesthetic changes: fixed spacings in csv table outputs
%
%4/27/16 PC
% added version number to table ouputs
% minor code changes to how cell array is created for datap when create
%  excel file for saving IDs
%
%7/18/16 PC
% fixed the legend of icc and betstddev plot
%
%7/19/16 PC
% removed 'n Inlcluded' from stddev stable to avoid confusion (as stddev
%  takes into account variance from entire sample)
% fixed naming of variable that wasn't correctly grabbing the number of
%  trials to plot for the depplot
%
%7/20/16 PC
% removed ICCs from overall dep table and put in stddev table
%  changes associated with viewing and saving the information
%
%7/23/16 PC
% updated commments (statistics and wavelet toolbox no longer required)
% pulled out the reliability calculations and put into their own function:
%  
%7/24/16 PC
% finished incoporating new dependability function (era_dep)
% incorporated new functions for displaying all figures and tables
% era_relfigures is not a wrapper for displaying (mostly meant to be used
%  by the gui)
%
%1/19/17 PC
% updated copyright
%8/17/17 PC
% cleaned up comments a bit

%somersault through inputs
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_relfigures for more information on optional inputs'));
    end
    
    %check whether era_data was provided
    %If it is not found, set display error.
    ind = find(strcmp('era_data',varargin),1);
    if ~isempty(ind)
        era_data = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: era_data not specified \n\n',... 
        'Please input the era_data to be loaded \n'));
    end
   
    %check if the dependability cutoff was specified 
    ind = find(strcmp('depcutoff',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            depcutoff = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            depcutoff = varargin{ind+1}; 
        end
    else 
        depcutoff = .70; %default level is .70
    end
    
    %check if plotdep is provided
    ind = find(strcmp('plotdep',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            pdep = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            pdep = varargin{ind+1}; 
        end
    else 
        pdep = 1; %default is 1
    end
    
    %check if picc is provided
    ind = find(strcmp('ploticc',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            picc = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            picc = varargin{ind+1}; 
        end
    else 
        picc = 1; %default is 1
    end

    %check if showinct is provided
    ind = find(strcmp('showinct',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            showinct = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            showinct = varargin{ind+1}; 
        end
    else 
        showinct = 1; %default is 1
    end
    
    %check if showoverallt is provided
    ind = find(strcmp('showoverallt',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            showoverallt = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            showoverallt = varargin{ind+1}; 
        end
    else 
        showoverallt = 1; %default is 1
    end
    
    %check if plotntrials is provided
    ind = find(strcmp('plotntrials',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotntrials = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            plotntrials = varargin{ind+1}; 
        end
    else 
        plotntrials = 50; %default is 50
    end
    
    %check if showstddevt is provided
    ind = find(strcmp('showstddevt',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            showstddevt = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            showstddevt = varargin{ind+1}; 
        end
    else 
        showstddevt = 1; %default is 1
    end
    
    %check if plotbetstddev is provided
    ind = find(strcmp('plotbetstddev',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotbetstddev = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            plotbetstddev = varargin{ind+1}; 
        end
    else 
        plotbetstddev = 1; %default is 1
    end
    
   	%check if plotwitstddev is provided
    ind = find(strcmp('plotwitstddev',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotwitstddev = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            plotwitstddev = varargin{ind+1}; 
        end
    else 
        plotwitstddev = 0; %default is 1
    end
    
    %check if plotdepline is provided
    ind = find(strcmp('plotdepline',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotdepline = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            plotdepline = varargin{ind+1}; 
        end
    else 
        plotdepline = 1; %default is 1 (lower limit of credible interval)
    end
    
    %check if meascutoff is provided
    ind = find(strcmp('meascutoff',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            meascutoff = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            meascutoff = varargin{ind+1}; 
        end
    else 
        meascutoff = 1; %default is 1 (lower limit of credible interval)
    end
    
    %check if depcentmeas is provided
    ind = find(strcmp('depcentmeas',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            depcentmeas = cell2mat(varargin{ind+1}); 
        elseif isnumeric(varargin{ind+1})
            depcentmeas = varargin{ind+1}; 
        end
    else 
        depcentmeas = 1; %default is 1 (mean)
    end
    
elseif ~isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help era_relfigures for more information on inputs'));
    
end %if ~isempty(varargin)

%calculate reliabitliy information to be used for plotting and tables
[era_data, relerr] = era_relsummary('era_data',era_data,'depcutoff',depcutoff,...
  'meascutoff',meascutoff,'depcentmeas',depcentmeas);

%if no good data were found then abort so you user can specify a different
%reliability threshold
if relerr.nogooddata == 1
    return;
end

%check which figures and tables have been specified to be shown
%prepare the figures and display them

%Plot that shows the relationship between the number of trials retained 
%for averaging and dependability
if pdep == 1
    era_depvtrialsplot('era_data',era_data,'trials',[1 plotntrials],...
        'depline',plotdepline,'depcutoff',depcutoff);
end

%Plot that compares the intraclass correlation coefficients for each
%group and/or condition
if picc == 1
    era_ptintervalplot('era_data',era_data,'stat','icc');
end

%plot that shows between-person standard deviation
if plotbetstddev == 1
    era_ptintervalplot('era_data',era_data,'stat','bet');
end

%plat that shows within-person standard deviation
if plotwitstddev == 1
    era_ptintervalplot('era_data',era_data,'stat','wit');
end

%table displaying trial cutoff information
if showinct == 1
    era_depcutofft('era_data',era_data,'gui',1);
end

%table displaying overall dependability information
if showoverallt == 1
    era_depoverallt('era_data',era_data,'gui',1);
end

%table including ICCs and between- and within-person standard deviations
if showstddevt == 1
    era_variancet('era_data',era_data,'gui',1);
end

%show an error if the dependability threshold was never met for the cutoff
if relerr.trlcutoff == 1
    errorstr = {};
    errorstr{end+1} = 'Trial cutoffs for adequate dependability could not be calculated';
    errorstr{end+1} = 'Data are too variable or there are not enough trials';

    errordlg(errorstr);
end

%show an error if none of the data had enough trials to meet the cutoff
if relerr.trlmax == 1
    errorstr = {};
    errorstr{end+1} = 'Not enough trials are present in the current data';
    errorstr{end+1} = 'Cutoffs represent an extrapolation beyond the data';
    
    errordlg(errorstr);
end

end

