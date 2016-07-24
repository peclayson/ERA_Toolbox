function era_relfigures(varargin)
%Creates various figures and tables for dependability data. See user manual
% for more specific information about each figure.
%
%era_relfigures('data',ERAData)
%
%Last Modified 7/24/16
%
%Required Inputs:
% data - structure array containing results of dependability analyses using
%  cmdstan. See era_computerel for more information.
%
%Optional Inputs:
% depcutoff - dependability level to use for cutoff when deciding the
%  minimum number of trials needed to achieve this specified level of
%  dependability
% plotdep - plot the table displaying dependability v number of trials
%  included in average
% ploticc - plot the intraclass correlation coefficients for data. This ICC
%  can be interpreted as the proportion of total variance that is between
%  persons.
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
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('data',varargin),1);
    if ~isempty(ind)
        REL = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
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
    
    %check if depplotntrials is provided
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

%create a data structure for storing outputs
data = struct;

%create variable to specify whether there are bad/unreliable data
%default state: 0, will be changed to 1 if there is a problem
poorrel = struct();
poorrel.trlcutoff = 0;
poorrel.trlmax = 0;

%check whether any groups exist
if strcmpi(REL.groups,'none')
    ngroups = 1;
    %gnames = cellstr(REL.groups);
    gnames ={''};
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

%check whether any events exist
if strcmpi(REL.events,'none')
    nevents = 1;
    %enames = cellstr(REL.events);
    enames = {''};
else
    nevents = length(REL.events);
    enames = REL.events(:);
end


%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

%extract information from REL and store in data for crunching
switch analysis
    case 1 %1 - no groups or event types to consider

        gloc = 1;
        eloc = 1;
        
        data.g(gloc).e(eloc).label = REL.out.labels(gloc);
        data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
        data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,gloc);
        data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,gloc);
        data.g(gloc).e(eloc).elabel = cellstr('none');
        data.g(gloc).glabel = gnames(gloc);
        
    case 2 %2 - possible multiple groups but no event types to consider 
        
        eloc = 1;
        
        for gloc=1:length(REL.out.labels)
            
            data.g(gloc).e(eloc).label = REL.out.labels(gloc);
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
            data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,gloc);
            data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,gloc);
            data.g(gloc).e(eloc).elabel = cellstr('none');
            data.g(gloc).glabel = gnames(gloc);
    
        end
        
    case 3 %3 - possible event types but no groups to consider
        
        gloc = 1;
        
        for eloc=1:length(REL.out.labels)

            data.g(gloc).e(eloc).label = REL.out.labels(eloc);
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
            data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,eloc);
            data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,eloc);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
    
        end
        
    case 4 %4 - possible groups and event types to consider
        for i=1:length(REL.out.labels)
                
            %use the underscores that were added in era_computerel to
            %differentiate where the group and event the data are for
            lblstr = strsplit(REL.out.labels{i},'_;_');

            eloc = find(ismember(enames,lblstr(1)));
            gloc = find(ismember(gnames,lblstr(2)));

%             if isempty(eloc); eloc = 1; end;
%             if isempty(gloc); gloc = 1; end;

            data.g(gloc).e(eloc).label = REL.out.labels(i);
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,i);
            data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,i);
            data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,i);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
    
        end
end %switch analysis

%store the cutoff in relsummary to pass to other functions more easily
relsummary.depcutoff = depcutoff;

%store which measure was used to specify cutoff
switch meascutoff
    case 1
        relsummary.meascutoff = 'Lower Limit of 95% Credible Interval';
    case 2
        relsummary.meascutoff = 'Point Estimate';
    case 3
        relsummary.meascutoff = 'Upper Limit of 95% Credible Interval';
end

%create an x-axis for the number of observations
ntrials = plotntrials;
x = 1:ntrials;

%create an empty array for storing information into
plotrel = zeros(ntrials,0);

%see how many subplots are needed
if nevents > 2
    xplots = ceil(sqrt(nevents));
    yplots = ceil(sqrt(nevents));
elseif nevents == 2
    xplots = 1;
    yplots = 2;
elseif nevents == 1
    xplots = 1;
    yplots = 1;
end

%create the figure
depplot = figure('Visible','Off');
set(gcf,'NumberTitle','Off');
depplot.Position = [125 630 900 450];
fsize = 16;

%extract the data and create the subplots for depplot
for eloc=1:nevents
    for gloc=1:ngroups
        switch plotdepline
            case 1 %lower limit
                [plotrel(x,gloc),~,~] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);
                plottitle = 'Lower Limit of 95% Credible Interval';
            case 2 %point estimate
                [~,plotrel(x,gloc),~] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);
                plottitle = 'Point Estimate';
            case 3 %upper limit
                [~,~,plotrel(x,gloc)] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);
                plottitle = 'Upper Limit of 95% Credible Interval';
        end
    end
    
    pref = 'Dependability v Number of Trials: ';
    set(gcf,'Name',[pref plottitle]);
    subplot(yplots,xplots,eloc); 
    h = plot(x,plotrel);
    
    axis([0 ntrials 0 1]);
    set(gca,'fontsize',16);
    
    if ~strcmpi(enames{eloc},'none')
        title(enames{eloc},'FontSize',20);
    end
    
    ylabel('Dependability','FontSize',fsize);
    xlabel('Number of Observations','FontSize',fsize);
    hline = refline(0,depcutoff);
    set(hline,'Color','b','LineStyle',':');
    if analysis ~= 1 && analysis ~= 3
        leg = legend(gnames{:},'Location','southeast');
        set(leg,'FontSize',fsize);
    end

end

%compute dependability data for each group and event
switch analysis
    case 1 %no groups or event types to consider
        
        %the same generic structure is used for relsummary, so the event
        %and group locations will both be 1 for the data
        eloc = 1;
        gloc = 1;
        
        %grab the group names (if present)
        relsummary.group(gloc).name = gnames{gloc};
        
        %set the event name as measure
        relsummary.group(gloc).event(eloc).name = 'measure';
        
        try
            %create empty arrays for storing dependability information
            trltable = varfun(@length,REL.data{1},'GroupingVariables',{'id'});
        catch
            %create empty arrays for storing dependability information
            trltable = varfun(@length,REL.data,'GroupingVariables',{'id'});
        end
        
        %compute dependabiltiy
        ntrials = max(trltable.GroupCount(:)) + 100;
        [llrel,mrel,ulrel] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);

        %find the number of trials to reach cutoff 
        switch meascutoff
            case 1 
                trlcutoff = find(llrel >= depcutoff, 1);
            case 2
                trlcutoff = find(mrel >= depcutoff, 1);
            case 3
                trlcutoff = find(ulrel >= depcutoff, 1);
        end
        
        %see whether the trial cutoff was found. If not store all values as
        %-1. If it is found, store the dependability information about the 
        %cutoffs.
        if isempty(trlcutoff)
            
            poorrel.trlcutoff = 1;
            trlcutoff = -1;
            relsummary.group(gloc).event(eloc).trlcutoff = -1;
            relsummary.group(gloc).event(eloc).rel.m = -1;
            relsummary.group(gloc).event(eloc).rel.ll = -1;
            relsummary.group(gloc).event(eloc).rel.ul = -1;
            
        elseif ~isempty(trlcutoff)
            
            %store information about cutoffs
            relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
            relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
            relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
            relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);

        end    
        
        %find the participants without enough trials based on the cutoffs
        if isempty(trlcutoff) %if the cutoff was not found in the data
            
            %Get information for the participants
            datatrls = REL.data;

            trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
            
            %define all of the data as bad, because the cutoff wasn't even
            %found
            relsummary.group(gloc).event(eloc).eventgoodids = 'none';
            relsummary.group(gloc).event(eloc).eventbadids =...
                trltable.id;

            relsummary.group(gloc).goodids = 'none';
            relsummary.group(gloc).badids = ...
                relsummary.group(gloc).event(eloc).eventbadids;
            
            %calculate the dependability estimates for the overall data
            datatable = REL.data;

            badids = table(relsummary.group(gloc).badids);

            trlcdata = innerjoin(datatable, badids,...
                'LeftKeys', 'id', 'RightKeys', 'Var1',...
                'LeftVariables', {'id' 'meas'});

            trltable = varfun(@length,trlcdata,...
                'GroupingVariables',{'id'});

            trlmean = mean(trltable.GroupCount);
            trlmed = median(trltable.GroupCount);

            %calculate dependability using either the mean or median (based
            %on user input)
            switch depcentmeas
                case 1
                    depcent = trlmean;
                case 2
                    depcent = trlmed;
            end
            
            
            [lldep,mdep,uldep] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',depcent,'CI',.95);

            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depcent;

            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            relsummary.group(gloc).event(eloc).goodn = 0;
            
        elseif ~isempty(trlcutoff)
            
            %Get trial information for those that meet cutoff
            datatrls = REL.data;

            trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

            ind2include = trltable.GroupCount >= trlcutoff;
            ind2exclude = trltable.GroupCount < trlcutoff;
            
            %store the good ids and the bad ids (ids that don't meet the
            %cutoff)
            relsummary.group(gloc).event(eloc).eventgoodids =...
                trltable.id(ind2include);
            relsummary.group(gloc).event(eloc).eventbadids =...
                trltable.id(ind2exclude);

            relsummary.group(gloc).goodids = ...
                relsummary.group(gloc).event(eloc).eventgoodids;
            relsummary.group(gloc).badids = ...
                relsummary.group(gloc).event(eloc).eventbadids;

            datatable = REL.data;
            
            goodids = table(relsummary.group(gloc).goodids);

            trlcdata = innerjoin(datatable, goodids,...
                'LeftKeys', 'id', 'RightKeys', 'Var1',...
                'LeftVariables', {'id' 'meas'});

            trltable = varfun(@length,trlcdata,...
                'GroupingVariables',{'id'});

            trlmean = mean(trltable.GroupCount);
            trlmed = median(trltable.GroupCount);

            %calculate dependability using either the mean or median (based
            %on user input)
            switch depcentmeas
                case 1
                    depcent = trlmean;
                case 2
                    depcent = trlmed;
            end

            [lldep,mdep,uldep] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',depcent,'CI',.95);

            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depcent;

            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            relsummary.group(gloc).event(eloc).goodn = height(goodids);
            
            %just in case a cutoff was extrapolated, the dependability
            %information will be overwritten
            if  relsummary.group(gloc).event(eloc).trlcutoff >...
                    relsummary.group(gloc).event(eloc).trlinfo.max
                
                poorrel.trlmax = 1;
                relsummary.group(gloc).event(eloc).rel.m = -1;
                relsummary.group(gloc).event(eloc).rel.ll = -1;
                relsummary.group(gloc).event(eloc).rel.ul = -1;

            end
            
            
        end
        
        %calculate iccs, between-subject standard deviations, and
        %within-subject standard deviations
        relsummary.group(gloc).event(eloc).icc.m = ...
            mean(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
            data.g(gloc).e(eloc).sig_e.raw)); 
        relsummary.group(gloc).event(eloc).icc.ll = ...
            quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
            data.g(gloc).e(eloc).sig_e.raw),.025); 
        relsummary.group(gloc).event(eloc).icc.ul = ...
            quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
            data.g(gloc).e(eloc).sig_e.raw),.975); 

        relsummary.group(gloc).event(eloc).betsd.m = ...
            mean(data.g(gloc).e(eloc).sig_u.raw); 
        relsummary.group(gloc).event(eloc).betsd.ll = ...
            quantile(data.g(gloc).e(eloc).sig_u.raw,.025); 
        relsummary.group(gloc).event(eloc).betsd.ul = ...
            quantile(data.g(gloc).e(eloc).sig_u.raw,.975);
        
        relsummary.group(gloc).event(eloc).witsd.m = ...
            mean(data.g(gloc).e(eloc).sig_e.raw); 
        relsummary.group(gloc).event(eloc).witsd.ll = ...
            quantile(data.g(gloc).e(eloc).sig_e.raw,.025); 
        relsummary.group(gloc).event(eloc).witsd.ul = ...
            quantile(data.g(gloc).e(eloc).sig_e.raw,.975);
       
        
    case 2 %possible multiple groups but no event types to consider
        
        %since the same generic structure is used for relsummary, the event
        %location will be defined as 1. 
       	eloc = 1;

        for gloc=1:ngroups %loop through each group

            %store the name of the group
            if eloc == 1
                relsummary.group(gloc).name = gnames{gloc};
            end

            %store the name of the event
            relsummary.group(gloc).event(eloc).name = enames{eloc};
            
            try
                %create empty arrays for storing dependability information
                trltable = varfun(@length,REL.data{1},...
                    'GroupingVariables',{'id'});
            catch
                %create empty arrays for storing dependability information
                trltable = varfun(@length,REL.data,...
                    'GroupingVariables',{'id'});
            end
        

            %compute dependability
            ntrials = max(trltable.GroupCount(:)) + 100;
            [llrel,mrel,ulrel] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);

            %find the number of trials to reach cutoff
            switch meascutoff
                case 1 
                    trlcutoff = find(llrel >= depcutoff, 1);
                case 2
                    trlcutoff = find(mrel >= depcutoff, 1);
                case 3
                    trlcutoff = find(ulrel >= depcutoff, 1);
            end
            
            %see whether the trial cutoff was found. If not store all 
            %values as -1. If it is found, store the dependability 
            %information about the cutoffs.
            
            if isempty(trlcutoff)

                poorrel.trlcutoff = 1;
                trlcutoff = -1;
                relsummary.group(gloc).event(eloc).trlcutoff = -1;
                relsummary.group(gloc).event(eloc).rel.m = -1;
                relsummary.group(gloc).event(eloc).rel.ll = -1;
                relsummary.group(gloc).event(eloc).rel.ul = -1;

            elseif ~isempty(trlcutoff)

                %store information about cutoffs
                relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);

            end 
            
            %if the trial cutoff was not found in the specified data
            if trlcutoff == -1
                
                %store all the ids as bad
                datatrls = REL.data;
                ind = strcmp(datatrls.group,gnames{gloc});
                datatrls = datatrls(ind,:);

                trltable = varfun(@length,datatrls,...
                    'GroupingVariables',{'id'});

                relsummary.group(gloc).event(eloc).eventgoodids = 'none';
                relsummary.group(gloc).event(eloc).eventbadids =...
                    trltable.id;   
                
            else
                
                %store the ids for participants with enough trials as good,
                %and the ids for participants with too few trials as bad
                datatrls = REL.data;
                ind = strcmp(datatrls.group,gnames{gloc});
                datatrls = datatrls(ind,:);

                trltable = varfun(@length,datatrls,...
                    'GroupingVariables',{'id'});

                ind2include = trltable.GroupCount >= trlcutoff;
                ind2exclude = trltable.GroupCount < trlcutoff;

                relsummary.group(gloc).event(eloc).eventgoodids =...
                    trltable.id(ind2include);
                relsummary.group(gloc).event(eloc).eventbadids =...
                    trltable.id(ind2exclude);
                
            end

            
        end

        %since there is only 1 event, those participants with good data for
        %the event will be store as participants having good data for the
        %whole group (this step makes more sense when there is more than
        %one event per group)
        
        for gloc=1:ngroups
        
            relsummary.group(gloc).goodids = ...
                relsummary.group(gloc).event(eloc).eventgoodids;
            relsummary.group(gloc).badids = ...
                relsummary.group(gloc).event(eloc).eventbadids;
            
        end
        
        %calculate the dependability for the overall data for each group
        for gloc=1:ngroups

            datatable = REL.data;
            ind = strcmp(datatable.group,gnames{gloc});
            datasubset = datatable(ind,:);
            
            %only factor in the trial counts from those with good data
            
            if ~strcmp(relsummary.group(gloc).goodids,'none')
                goodids = table(relsummary.group(gloc).goodids);
            elseif strcmp(relsummary.group(gloc).goodids,'none')
                goodids = table(relsummary.group(gloc).goodids);
            end

            trlcdata = innerjoin(datasubset, goodids,...
                'LeftKeys', 'id', 'RightKeys', 'Var1',...
                'LeftVariables', {'id' 'meas'});

            trltable = varfun(@length,trlcdata,...
                'GroupingVariables',{'id'});

            trlmean = mean(trltable.GroupCount);
            trlmed = median(trltable.GroupCount);

            %calculate dependability using either the mean or median (based
            %on user input)
            switch depcentmeas
                case 1
                    depcent = trlmean;
                case 2
                    depcent = trlmed;
            end
            
            [lldep,mdep,uldep] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',depcent,'CI',.95);

            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depcent;
            
            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            
            %just in case a cutoff was extrapolated, the dependability
            %information will be overwritten
            if  relsummary.group(gloc).event(eloc).trlcutoff >...
                    relsummary.group(gloc).event(eloc).trlinfo.max
                
                poorrel.trlmax = 1;
                relsummary.group(gloc).event(eloc).rel.m = -1;
                relsummary.group(gloc).event(eloc).rel.ll = -1;
                relsummary.group(gloc).event(eloc).rel.ul = -1;

            end
            
            if ~strcmp(relsummary.group(gloc).goodids,'none')
                relsummary.group(gloc).event(eloc).goodn = height(goodids);
            elseif strcmp(relsummary.group(gloc).goodids,'none')
                relsummary.group(gloc).event(eloc).goodn = 0;
            end
            
            %calculate iccs, between-subject standard deviations, and
            %within-subject standard deviations            
            
            relsummary.group(gloc).event(eloc).icc.m = ...
                mean(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw)); 
            relsummary.group(gloc).event(eloc).icc.ll = ...
                quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw),.025); 
            relsummary.group(gloc).event(eloc).icc.ul = ...
                quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw),.975);
            
            relsummary.group(gloc).event(eloc).betsd.m = ...
                mean(data.g(gloc).e(eloc).sig_u.raw); 
            relsummary.group(gloc).event(eloc).betsd.ll = ...
                quantile(data.g(gloc).e(eloc).sig_u.raw,.025); 
            relsummary.group(gloc).event(eloc).betsd.ul = ...
                quantile(data.g(gloc).e(eloc).sig_u.raw,.975);

            relsummary.group(gloc).event(eloc).witsd.m = ...
                mean(data.g(gloc).e(eloc).sig_e.raw); 
            relsummary.group(gloc).event(eloc).witsd.ll = ...
                quantile(data.g(gloc).e(eloc).sig_e.raw,.025); 
            relsummary.group(gloc).event(eloc).witsd.ul = ...
                quantile(data.g(gloc).e(eloc).sig_e.raw,.975);
       
            
        end
        
        
    case 3 %possible event types but no groups to consider
        
        %since the relsummary structure is generic for any number of groups
        %or events, the data will count as being from 1 group. 
        gloc = 1;
        
        %cylce through each event
        for eloc=1:nevents
                
            %get the group name
            if eloc == 1
                relsummary.group(gloc).name = gnames{gloc};
            end
            
            %get the event name
            relsummary.group(gloc).event(eloc).name = enames{eloc};
            
            %create empty arrays for storing dependability information
            trltable = varfun(@length,REL.data{1},...
                'GroupingVariables',{'id'});
            
            %compute dependability
            ntrials = max(trltable.GroupCount(:)) + 100;
            [llrel,mrel,ulrel] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);

            %find the number of trials to reach cutoffl
            switch meascutoff
                case 1 
                    trlcutoff = find(llrel >= depcutoff, 1);
                case 2
                    trlcutoff = find(mrel >= depcutoff, 1);
                case 3
                    trlcutoff = find(ulrel >= depcutoff, 1);
            end

            
            if isempty(trlcutoff) %if a cutoff wasn't found

                poorrel.trlcutoff = 1;
                trlcutoff = -1;
                relsummary.group(gloc).event(eloc).trlcutoff = -1;
                relsummary.group(gloc).event(eloc).rel.m = -1;
                relsummary.group(gloc).event(eloc).rel.ll = -1;
                relsummary.group(gloc).event(eloc).rel.ul = -1;

                datatrls = REL.data{eloc};

                trltable = varfun(@length,datatrls,...
                    'GroupingVariables',{'id'});

                %store all of the participant ids as bad, becuase none
                %reached the cutoff
                relsummary.group(gloc).event(eloc).eventgoodids =...
                    'none';
                relsummary.group(gloc).event(eloc).eventbadids =...
                    trltable.id;

            elseif ~isempty(trlcutoff) %if a cutoff was found

                %store information about cutoffs
                relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);

                datatrls = REL.data{eloc};

                trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                ind2include = trltable.GroupCount >= trlcutoff;
                ind2exclude = trltable.GroupCount < trlcutoff;
        
                %store which participants had enough trials to meet cutoff
                relsummary.group(gloc).event(eloc).eventgoodids =...
                    trltable.id(ind2include);
                relsummary.group(gloc).event(eloc).eventbadids =...
                    trltable.id(ind2exclude);

            end 
                
        end
        
        %cycle through the events and store information about which
        %participants have good data across all event types.
        tempids = {};
        badids = [];
        for eloc=1:nevents
            if ~strcmp(relsummary.group(gloc).event(eloc).eventgoodids,'none')
                tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                if eloc == 1
                    badids = relsummary.group(gloc).event(eloc).eventbadids;
                elseif eloc > 1
                    [~,ind]=setdiff(tempids{1},tempids{end});
                    new = tempids{1};
                    badids = vertcat(badids,...
                            relsummary.group(gloc).event(eloc).eventbadids,...
                            new(ind));
                    new(ind) = [];
                    tempids{1} = new;
                end
            end
        end
        
        %if none the participants had good data for all events
        if isempty(tempids)
            tempids{1} = 'none';
        end
        
        %store information about participants with good/bad ids
        relsummary.group(gloc).goodids = tempids{1};
        relsummary.group(gloc).badids = unique(badids);
        
        %cycle through each event
        for eloc=1:nevents
                
            %get trial information to use for dependability estimates. only
            %participants with enough data will contribute toward trial
            %counts
            datatable = REL.data{eloc};

            if ~strcmp(relsummary.group(gloc).goodids,'none')
                goodids = table(relsummary.group(gloc).goodids);
            else
                goodids = table(relsummary.group(gloc).badids);
            end

            trlcdata = innerjoin(datatable, goodids,...
                'LeftKeys', 'id', 'RightKeys', 'Var1',...
                'LeftVariables', {'id' 'meas'});

            trltable = varfun(@length,trlcdata,...
                'GroupingVariables',{'id'});

            trlmean = mean(trltable.GroupCount);
            trlmed = median(trltable.GroupCount);

            %calculate dependability using either the mean or median (based
            %on user input)
            switch depcentmeas
                case 1
                    depcent = trlmean;
                case 2
                    depcent = trlmed;
            end

            [lldep,mdep,uldep] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',depcent,'CI',.95);
                
            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depcent;

            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            
            %just in case a cutoff was extrapolated, the dependability
            %information will be overwritten
            if  relsummary.group(gloc).event(eloc).trlcutoff >...
                    relsummary.group(gloc).event(eloc).trlinfo.max

                poorrel.trlmax = 1;
                relsummary.group(gloc).event(eloc).rel.m = -1;
                relsummary.group(gloc).event(eloc).rel.ll = -1;
                relsummary.group(gloc).event(eloc).rel.ul = -1;

            end

            if ~strcmp(relsummary.group(gloc).goodids,'none')
                relsummary.group(gloc).event(eloc).goodn = height(goodids);
            else
                relsummary.group(gloc).event(eloc).goodn = 0;
            end

            %calculate iccs, between-subject standard deviations, and
            %within-subject standard deviations
            relsummary.group(gloc).event(eloc).icc.m = ...
                mean(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw)); 
            relsummary.group(gloc).event(eloc).icc.ll = ...
                quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw),.025); 
            relsummary.group(gloc).event(eloc).icc.ul = ...
                quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw),.975);

            relsummary.group(gloc).event(eloc).betsd.m = ...
                mean(data.g(gloc).e(eloc).sig_u.raw); 
            relsummary.group(gloc).event(eloc).betsd.ll = ...
                quantile(data.g(gloc).e(eloc).sig_u.raw,.025); 
            relsummary.group(gloc).event(eloc).betsd.ul = ...
                quantile(data.g(gloc).e(eloc).sig_u.raw,.975);

            relsummary.group(gloc).event(eloc).witsd.m = ...
                mean(data.g(gloc).e(eloc).sig_e.raw); 
            relsummary.group(gloc).event(eloc).witsd.ll = ...
                quantile(data.g(gloc).e(eloc).sig_e.raw,.025); 
            relsummary.group(gloc).event(eloc).witsd.ul = ...
                quantile(data.g(gloc).e(eloc).sig_e.raw,.975);
                
        end
                
    case 4 %groups and event types to consider
        
        %cycle through each event
        for eloc=1:nevents
            
            %cycle through each group
            for gloc=1:ngroups
                
                %store the group name
                if eloc == 1
                    relsummary.group(gloc).name = gnames{gloc};
                end
               
                %store the event name
                relsummary.group(gloc).event(eloc).name = enames{eloc};
                
                %create empty arrays for storing dependability information
                trltable = varfun(@length,REL.data{1},...
                    'GroupingVariables',{'id'});

                %compute dependability
                ntrials = max(trltable.GroupCount(:)) + 100;
                [llrel,mrel,ulrel] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',[1 ntrials],'CI',.95);
                
                %find the number of trials to reach cutoff
                switch meascutoff
                    case 1 
                        trlcutoff = find(llrel >= depcutoff, 1);
                    case 2
                        trlcutoff = find(mrel >= depcutoff, 1);
                    case 3
                        trlcutoff = find(ulrel >= depcutoff, 1);
                end
                
                %if a cutoff was not found
                if isempty(trlcutoff) 
                
                    %store all the participant ids as having bad data
                    poorrel.trlcutoff = 1;
                    trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).rel.m = -1;
                    relsummary.group(gloc).event(eloc).rel.ll = -1;
                    relsummary.group(gloc).event(eloc).rel.ul = -1;
                    
                    datatrls = REL.data{eloc};
                    ind = strcmp(datatrls.group,gnames{gloc});
                    datatrls = datatrls(ind,:);
                
                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        'none';
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id;
                    
                elseif ~isempty(trlcutoff) %if a cutoff was found

                    %store information about cutoffs
                    relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                    relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);

                    datatrls = REL.data{eloc};
                    ind = strcmp(datatrls.group,gnames{gloc});
                    datatrls = datatrls(ind,:);
                    
                    %find ids with enough trials based on cutoff
                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                    ind2include = trltable.GroupCount >= trlcutoff;
                    ind2exclude = trltable.GroupCount < trlcutoff;

                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        trltable.id(ind2include);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id(ind2exclude);
                    
                end 

            end
        end
        
        %find the ids that have enough trials for each event type
        
        for gloc=1:ngroups 
            tempids = {};
            badids = [];
            
            for eloc=1:nevents
                
                if ~strcmp(relsummary.group(gloc).event(eloc).eventgoodids,'none')
                    tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                    
                    if eloc == 1
                        badids = relsummary.group(gloc).event(eloc).eventbadids;
                    elseif eloc > 1
                        [~,ind]=setdiff(tempids{1},tempids{end});
                        new = tempids{1};
                        badids = vertcat(badids,...
                            relsummary.group(gloc).event(eloc).eventbadids,...
                            new(ind));
                        new(ind) = [];
                        tempids{1} = new;
                    end
                    
                end
            end
            
            %if there were no participants with enough trials for all event
            %types within a group, store 'none'
            if isempty(tempids)
                tempids{1} = 'none';
            end
            
            relsummary.group(gloc).goodids = tempids{1};
            relsummary.group(gloc).badids = unique(badids);
            
        end
         
        %calculate dependability estimates for the overall data for each
        %group and event type
        for eloc=1:nevents

            for gloc=1:ngroups
                
                %pull the data to calculate the trial information from
                %participants with good data
                datatable = REL.data{eloc};
                ind = strcmp(datatable.group,gnames{gloc});
                datasubset = datatable(ind,:);
                
                if ~strcmp(relsummary.group(gloc).goodids,'none')
                    goodids = table(relsummary.group(gloc).goodids);
                else
                    goodids = table(relsummary.group(gloc).badids);
                end
                
                trlcdata = innerjoin(datasubset, goodids,...
                    'LeftKeys', 'id', 'RightKeys', 'Var1',...
                    'LeftVariables', {'id' 'meas'});
                
                trltable = varfun(@length,trlcdata,...
                    'GroupingVariables',{'id'});
                
                trlmean = mean(trltable.GroupCount);
                trlmed = median(trltable.GroupCount);

                %calculate dependability using either the mean or median 
                %(based on user input)
                switch depcentmeas
                    case 1
                        depcent = trlmean;
                    case 2
                        depcent = trlmed;
                end

                [lldep,mdep,uldep] = era_dep(...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,...
                    'obs',depcent,'CI',.95);
                
                relsummary.group(gloc).event(eloc).dep.m = mdep;
                relsummary.group(gloc).event(eloc).dep.ll = lldep;
                relsummary.group(gloc).event(eloc).dep.ul = uldep;
                relsummary.group(gloc).event(eloc).dep.meas = depcent;
                
                relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                
                %check if trial cutoff for an event exceeded the total
                %number of trials for any subjects
                if  relsummary.group(gloc).event(eloc).trlcutoff >...
                        relsummary.group(gloc).event(eloc).trlinfo.max
                    
                    poorrel.trlmax = 1;
                    relsummary.group(gloc).event(eloc).rel.m = -1;
                    relsummary.group(gloc).event(eloc).rel.ll = -1;
                    relsummary.group(gloc).event(eloc).rel.ul = -1;
                
                end

                if ~strcmp(relsummary.group(gloc).event(eloc).eventgoodids,'none')
                    relsummary.group(gloc).event(eloc).goodn = height(goodids);
                else
                    relsummary.group(gloc).event(eloc).goodn = 0;
                end
                
                %calculate iccs, between-subject standard deviations, and
                %within-subject standard deviations
                relsummary.group(gloc).event(eloc).icc.m = ...
                    mean(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw)); 
                relsummary.group(gloc).event(eloc).icc.ll = ...
                    quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw),.025); 
                relsummary.group(gloc).event(eloc).icc.ul = ...
                    quantile(era_icc(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw),.975); 
               
                relsummary.group(gloc).event(eloc).betsd.m = ...
                    mean(data.g(gloc).e(eloc).sig_u.raw); 
                relsummary.group(gloc).event(eloc).betsd.ll = ...
                    quantile(data.g(gloc).e(eloc).sig_u.raw,.025); 
                relsummary.group(gloc).event(eloc).betsd.ul = ...
                    quantile(data.g(gloc).e(eloc).sig_u.raw,.975);

                relsummary.group(gloc).event(eloc).witsd.m = ...
                    mean(data.g(gloc).e(eloc).sig_e.raw); 
                relsummary.group(gloc).event(eloc).witsd.ll = ...
                    quantile(data.g(gloc).e(eloc).sig_e.raw,.025); 
                relsummary.group(gloc).event(eloc).witsd.ul = ...
                    quantile(data.g(gloc).e(eloc).sig_e.raw,.975);
       
            end
        end
                
        
end %switch analysis

%store relsummary information
RELout = REL;
RELout.relsummary = relsummary;

%if the user wants an icc plot
if picc == 1
    iccplot = figure;
    iccplot.Position = [125 65 900 450];
    set(gcf,'NumberTitle','Off');
    set(gcf,'Name','ICC Estimates');
    miccm = zeros(nevents,ngroups);
    lliccm = zeros(nevents,ngroups);
    uliccm = zeros(nevents,ngroups);
    offsetm = zeros(nevents,ngroups);
    
    %grab icc information for each event and group
    for gloc=1:ngroups
       for eloc=1:nevents
           miccm(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.m;
           lliccm(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.m...
               - relsummary.group(gloc).event(eloc).icc.ll;
           uliccm(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.ul...
               - relsummary.group(gloc).event(eloc).icc.m;
           
           %figure out spacing for plot
           if gloc < median(1:ngroups)
               offsetm(eloc,gloc) = eloc - (.4/ngroups);
           elseif gloc > median(1:ngroups)
               offsetm(eloc,gloc) = eloc + (.4/ngroups);
           elseif gloc == median(1:ngroups)
               offsetm(eloc,gloc) = eloc;
           end
           
       end
    end
    
    %find the dimensions
    [e,g] = size(miccm); 
    
    %plot
    if ~(g > 1 && e == 1)
        iccplot = errorbar(offsetm,miccm,lliccm,uliccm,'Marker','.',...
            'MarkerSize',15,'LineWidth',1);
    elseif g > 1 && e == 1
        hold on
        for i = 1:g
            iccplot = errorbar(offsetm(i),miccm(i),lliccm(i),uliccm(i),...
                'Marker','.','MarkerSize',15,'LineWidth',1,...
                'DisplayName',gnames{i});
            legend('-DynamicLegend');
        end
        hold off
    end
    
    %add y-axis label
    ylabel('Intraclass Correlation Coefficient','FontSize',fsize);
    
    %remove extra lines
    if nevents > 1
        for i = 1:length(iccplot) 
            iccplot(i).LineStyle = 'none';
        end
    end
    
    %fix axes
    iccplot(1).Parent.XLim = [0 nevents+1+.25];
    if nevents == 1
        iccplot(1).Parent.YLim = [0 max(uliccm+miccm)+.05];
    elseif nevents > 1
        iccplot(1).Parent.YLim = [0 max(max(uliccm+miccm))+.05];
        xlabel('Event','FontSize',fsize);
    end
    
    iccplot(1).Parent.XTick = 1:nevents;
    iccplot(1).Parent.XTickLabel = enames;
    iccplot(1).Parent.FontSize = fsize;
    
    %add names to legend
    if ~(g > 1 && e == 1)
        pl = legend(iccplot);
        pl.String = gnames;
    end
    
    %change axis location and rotate plot
    iccplot(1).Parent.YAxisLocation = 'right';
    camroll(-90);
    
end

%if the user wants a between-person standard deviation plot
if plotbetstddev == 1
    sdplot = figure;
    sdplot.Position = [225 165 900 450];
    set(gcf,'NumberTitle','Off');
    set(gcf,'Name','Between-Person Standard Deviations');
    
    %create placeholders
    msd = zeros(nevents,ngroups);
    llsd = zeros(nevents,ngroups);
    ulsd = zeros(nevents,ngroups);
    offsetm = zeros(nevents,ngroups);
    
    %store between-subject std dev info
    for gloc=1:ngroups
       for eloc=1:nevents
           msd(eloc,gloc) = relsummary.group(gloc).event(eloc).betsd.m;
           llsd(eloc,gloc) = relsummary.group(gloc).event(eloc).betsd.m...
               - relsummary.group(gloc).event(eloc).betsd.ll;
           ulsd(eloc,gloc) = relsummary.group(gloc).event(eloc).betsd.ul...
               - relsummary.group(gloc).event(eloc).betsd.m;
           
           %figure out spacing
           if gloc < median(1:ngroups)
               offsetm(eloc,gloc) = eloc - (.4/ngroups);
           elseif gloc > median(1:ngroups)
               offsetm(eloc,gloc) = eloc + (.4/ngroups);
           elseif gloc == median(1:ngroups)
               offsetm(eloc,gloc) = eloc;
           end
           
       end
    end
    
    %find dimensions
    [e,g] = size(miccm); 
    
    %plot std devs
    if ~(g > 1 && e == 1)
        sdplot = errorbar(offsetm,msd,llsd,ulsd,'Marker','.',...
            'MarkerSize',15,'LineWidth',1);
    elseif g > 1 && e == 1
        hold on
        for i = 1:g
            sdplot = errorbar(offsetm(i),msd(i),llsd(i),ulsd(i),...
                'Marker','.','MarkerSize',15,'LineWidth',1,...
                'DisplayName',gnames{i});
            legend('-DynamicLegend');
        end
        hold off
    end

    %add y-axis label
    ylabel('Between-Person Standard Deviation','FontSize',fsize);
    
    %remove extra lines
    if nevents > 1
        for i = 1:length(sdplot)
            sdplot(i).LineStyle = 'none';
        end
    end
    
    %fix axes
    sdplot(1).Parent.XLim = [0 nevents+1+.25];
    if nevents == 1
        sdplot(1).Parent.YLim = [min(msd-llsd)-1.5 max(ulsd+msd)+1.5];
    elseif nevents > 1
        sdplot(1).Parent.YLim = [min(min(msd-llsd))-1.5 max(max(ulsd+msd))+1.5];
        xlabel('Event','FontSize',fsize);
    end
    
    sdplot(1).Parent.XTick = 1:nevents;
    sdplot(1).Parent.XTickLabel = enames;
    sdplot(1).Parent.FontSize = fsize;
    
    %add names to legend
    if ~(g > 1 && e == 1)
        pl = legend(sdplot);
        pl.String = gnames;
    end
    
    %change axis location and rotate plot
    sdplot(1).Parent.YAxisLocation = 'right';
    camroll(-90);
    
end

%reminder....
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

%create placeholders for displaying data in tables in guis
label = {};
trlcutoff = {};
dep = {};
overalldep = {};
mintrl = {};
maxtrl = {};
meantrl = {};
medtrl = {};
goodn = {};
badn = {};
icc = {};
betsd = {};
witsd = {};

%put data together to display in tables
for gloc=1:ngroups
    for eloc=1:nevents
        
        %label for group and/or event
        switch analysis
            case 1
                label{end+1} = 'Measurement';
            case 2
                label{end+1} = gnames{gloc};
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end
        
        %pull the trial cutoff
        trlcutoff{end+1} = relsummary.group(gloc).event(eloc).trlcutoff;
        
        %create a string with the dependability point estimate and credible
        %interval for cutoff data
        dep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).rel.m,...
            relsummary.group(gloc).event(eloc).rel.ll,...
            relsummary.group(gloc).event(eloc).rel.ul);
        
        %create a string with the dependability point estimate and credible
        %interval for overall data
        overalldep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).dep.m,...
            relsummary.group(gloc).event(eloc).dep.ll,...
            relsummary.group(gloc).event(eloc).dep.ul);
        
        %put together trial summary information
        mintrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.min;
        maxtrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.max;
        meantrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.mean;
        medtrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.med;
        
        %pull good and bad ns
        goodn{end+1} = relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = length(relsummary.group(gloc).badids);
        
        %create a string with the icc point estimate and credible interval
        icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).icc.m,...
            relsummary.group(gloc).event(eloc).icc.ll,...
            relsummary.group(gloc).event(eloc).icc.ul);
        
        %create a string with the between-person standard devation point 
        %estimate and credible interval
        betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).betsd.m,...
            relsummary.group(gloc).event(eloc).betsd.ll,...
            relsummary.group(gloc).event(eloc).betsd.ul);
        
        %create a string with the within-person standard devation point 
        %estimate and credible interval
        witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).witsd.m,...
            relsummary.group(gloc).event(eloc).witsd.ll,...
            relsummary.group(gloc).event(eloc).witsd.ul);
        
    end 
end

%create table for displaying between-person and within-person standard
%deviation information
stddevtable = table(label',betsd',witsd',icc');

stddevtable.Properties.VariableNames = {'Label',...
    'Between_StdDev','Within_StdDev','ICC'};

%create table to display the cutoff information with its associated
%dependability info
inctrltable = table(label',trlcutoff',dep');

inctrltable.Properties.VariableNames = {'Label','Trial_Cutoff',...
    'Dependability'};

%create table to describe the data including all trials 
overalltable = table(label',goodn',badn',overalldep',meantrl',...
    medtrl',mintrl',maxtrl');

overalltable.Properties.VariableNames = {'Label', ...
    'n_Included','n_Excluded', ...
    'Dependability', 'Mean_Num_Trials', 'Med_Num_Trials',...
    'Min_Num_Trials',...
    'Max_Num_Trials'};



%define parameters for figure size
figwidth = 674;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;
name = ['Point and 95% Interval Estimates for the Between-'...
    'and Within-Person Standard Deviations and ICCs'];

%create gui for standard-deviation table
era_stddev= figure('unit','pix','Visible','off',...
  'position',[1250 600 figwidth figheight],...
  'menub','no',...
  'name',name,...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(era_stddev,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    'Between- and Within-Person Standard Deviations and ICCs',...
    'Position',[0 row figwidth 25]);          

%Start a table
t = uitable('Parent',era_stddev,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(stddevtable));
set(t,'ColumnName',{'Label' 'Between Std Dev'...
    'Within Std Dev' 'ICC'});
set(t,'ColumnWidth',{200 140 140 140});
set(t,'RowName',[]);

%Create a save button that will take save the table
uicontrol(era_stddev,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/8 25 figwidth/4 50],...
    'Callback',{@era_savestddevtable,RELout,stddevtable}); 




%define parameters for figure size
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*1.6;

%create a gui to display the cutoff information table
era_inctrl= figure('unit','pix','Visible','off',...
  'position',[1150 700 figwidth figheight],...
  'menub','no',...
  'name',...
  'Results of Increasing Trials the Number of Trials on Dependability',...
  'numbertitle','off',...
  'resize','off');

%Add table title
uicontrol(era_inctrl,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    sprintf('Dependability Analyses, %0.2f Cutoff\n',...
    RELout.relsummary.depcutoff),...
    'Position',[0 row figwidth 20]);   

uicontrol(era_inctrl,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    sprintf('Cutoff Used the %s',...
    RELout.relsummary.meascutoff),...
    'Position',[0 row-20 figwidth 20]); 

%Start a table
t = uitable('Parent',era_inctrl,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(inctrltable));
set(t,'ColumnName',{'Label' 'Trial Cutoff' 'Dependability'});
set(t,'ColumnWidth',{200 'auto' 170});
set(t,'RowName',[]);

%Create a save button that will take save the table
uicontrol(era_inctrl,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/8 25 figwidth/4 50],...
    'Callback',{@era_saveinctrltable,RELout,inctrltable}); 



%define parameters for figure size
figwidth = 725;
figheight = 500;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%create a gui for displaying the overall trial information
era_overall= figure('unit','pix','Visible','off',...
  'position',[1150 150 figwidth figheight],...
  'menub','no',...
  'name','Dependability Analyses Including All Trials',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(era_overall,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Overall Dependability',...
    'Position',[0 row figwidth 25]);          

%Start a table
t = uitable('Parent',era_overall,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(overalltable));
set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
    'Dependability' 'Mean # of Trials' 'Med # of Trials'...
    'Min # of Trials' 'Max # of Trials'});
set(t,'ColumnWidth',{'auto' 'auto' 'auto' 110 'auto' 'auto' 'auto'});
set(t,'RowName',[]);

%Create a save button that will take save the table
uicontrol(era_overall,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/8 25 figwidth/4 50],...
    'Callback',{@era_saveoveralltable,RELout,overalltable}); 

%Create button that will save good/bad ids
uicontrol(era_overall,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save IDs',...
    'Position', [5*figwidth/8 25 figwidth/4 50],...
    'Callback',{@era_saveids,RELout}); 

%show the standard deviation table if it is wanted
if showstddevt == 1
    set(era_stddev,'Visible','on');
end

%show the dependability plot if it is wanted
if pdep == 1
    set(depplot,'Visible','on');
end

%show the cutoff table if it is wanted
if showinct == 1
    set(era_inctrl,'Visible','on');
end

%show the overall table if it is wanted
if showoverallt == 1
    set(era_overall,'Visible','on');
end

%show an error if the dependability threshold was never met for the cutoff
if poorrel.trlcutoff == 1
    errorstr = {};
    errorstr{end+1} = 'Trial cutoffs for adequate dependability could not be calculated';
    errorstr{end+1} = 'Data are too variable or there are not enough trials';

    errordlg(errorstr);
end

%show an error if none of the data had enough trials to meet the cutoff
if poorrel.trlmax == 1
    errorstr = {};
    errorstr{end+1} = 'Not enough trials are present in the current data';
    errorstr{end+1} = 'Cutoffs represent an extrapolation beyond the data';
    
    errordlg(errorstr);
end

end

function era_saveoveralltable(varargin)
%if the button to save the overall trial information table was pressed

%parse inputs
REL = varargin{3};
overalltable = varargin{4};

%ask the user where the file should be saved
if ~ismac %macs can't use xlswrite
    [savename, savepath] = uiputfile(...
        {'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
else
    [savename, savepath] = uiputfile(...
        {'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

%save either an excel or csv file
if strcmp(ext,'.xlsx')
    
    %print header information about the dataset
    filehead = {'Dependability Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('ERA Toolbox v%s',REL.eraver);
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',REL.filename);
    filehead{end+1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    filehead{end+1} = sprintf('Cutoff Threshold used the %s',...
        REL.relsummary.meascutoff);
    filehead{end+1} = sprintf('Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    filehead{end+1}='';
    filehead{end+1}='';
    
    %write table
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(overalltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    %print header information about dataset
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s','Dependability Table Generated on ');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',REL.eraver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'Cutoff Threshold used the %s\n',...
        REL.relsummary.meascutoff);
    fprintf(fid, 'Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s', strcat('Label,N Included,N Excluded,',...
        'Dependability,Mean Num of Trials,Med Num of Trials,',...
        'Min Num of Trials,Max Num of Trials'));
    fprintf(fid,' \n');
    
    %write the table information
    for i = 1:height(overalltable)
         formatspec = '%s,%d,%d,%s,%0.2f,%d,%d,%d\n';
         fprintf(fid,formatspec,char(overalltable{i,1}),...
             cell2mat(overalltable{i,2}),cell2mat(overalltable{i,3}),...
             char(overalltable{i,4}),cell2mat(overalltable{i,5}),...
             cell2mat(overalltable{i,6}),cell2mat(overalltable{i,7}),...
             cell2mat(overalltable{i,8}));
    end
    
    fclose(fid);
    
end

end

function era_saveids(varargin)
%if the user pressed the button to save which ids were considered good and
%which were considered bad (based on whether data met the cutoff thresholds

REL = varargin{3};

%ask the user where the file should be saved
if ~ismac
    [savename, savepath] = uiputfile(...
        {'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save the data?');
else
    [savename, savepath] = uiputfile(...
        {'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save the data?');
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

%save the information in either an excel or csv format
if strcmp(ext,'.xlsx')
    
    datap{1,1} = 'Data Generated on';
    datap{end+1,1} = datestr(clock); 
    datap{end+1,1} = sprintf('ERA Toolbox v%s',REL.eraver);
    datap{end+1,1} = '';
    datap{end+1,1} = sprintf('Dataset: %s',REL.filename);
    datap{end+1,1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    datap{end+1,1} = sprintf('Cutoff Threshold used the %s',...
        REL.relsummary.meascutoff);
    datap{end+1,1} = sprintf('Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    datap{end+1,1}='';
    datap{end+1,1}='';
    datap{end+1,1} = 'Good IDs'; datap{end,2} = 'Bad IDs';
    
    gids = [];
    bids = [];
    srow = length(datap);
    
    for j=1:length(REL.relsummary.group)
        gids = [gids;REL.relsummary.group(j).goodids(:)];
        bids = [bids;REL.relsummary.group(j).badids(:)];
    end
    
    for i = 1:length(gids)
        datap{i+srow,1}=char(gids(i));
    end

    for i = 1:length(bids)
        datap{i+srow,2}=char(bids(i));
    end

    xlswrite(fullfile(savepath,savename),datap);
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Data Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',REL.eraver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'Cutoff Threshold used the %s\n',...
        REL.relsummary.meascutoff);
    fprintf(fid,'Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    fprintf(fid,'%s\n','Good IDs,Bad IDs');
    
    gids = [];
    bids = [];
    
    for j=1:length(REL.relsummary.group)
        gids = [gids;REL.relsummary.group(j).goodids(:)];
        bids = [bids;REL.relsummary.group(j).badids(:)];
    end

    maxlength = max([length(gids) length(bids)]);
    minlength = min([length(gids) length(bids)]);
    
    if minlength == length(gids)
        whichlonger = 1;
    elseif maxlength == length(gids)
        whichlonger = 2;
    end
    
    for i = 1:maxlength
        if i <= minlength
            fprintf(fid,'%s,%s\n',gids{i},bids{i});
        elseif i > minlength && whichlonger == 1
            fprintf(fid,',%s\n',bids{i});
        elseif i > minlength && whichlonger == 2
            fprintf(fid,'%s,\n',gids{i});
        end
    end
    
    fclose(fid);
    
end


end

function era_saveinctrltable(varargin)
%if the user pressed the button to save the table with the information
%about the cutoffs

%parse inputs
REL = varargin{3};
incltrltable = varargin{4};

%ask the user where the file should be saved
if ~ismac %macs can't use xlswrite
    [savename, savepath] = uiputfile(...
        {'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
else
    [savename, savepath] = uiputfile(...
        {'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

%save as either excel or csv file
if strcmp(ext,'.xlsx')
    
    filehead = {'Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('ERA Toolbox v%s',REL.eraver);
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',REL.filename);
    filehead{end+1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    filehead{end+1} = sprintf('Cutoff Threshold used the %s',...
        REL.relsummary.meascutoff);
    filehead{end+1} = sprintf('Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(incltrltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',REL.eraver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'Cutoff Threshold used the %s\n',...
        REL.relsummary.meascutoff);
    fprintf(fid,'Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s', strcat('Label,Trial Cutoff,Dependability'));
    fprintf(fid,' \n');
    
    for i = 1:height(incltrltable)
         formatspec = '%s,%d,%s\n';
         fprintf(fid,formatspec,char(incltrltable{i,1}),...
             cell2mat(incltrltable{i,2}),char(incltrltable{i,3}));
    end
    
    fclose(fid);
    
end


end


function era_savestddevtable(varargin)
%if the user pressed the button to save the table with the information
%about between- and within-person standard deviations

%parse inputs
REL = varargin{3};
stddevtable = varargin{4};

%ask the user where the file should be saved
if ~ismac %macs can't use xlswrite
    [savename, savepath] = uiputfile(...
        {'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
else
    [savename, savepath] = uiputfile(...
        {'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

%save as either excel or csv file
if strcmp(ext,'.xlsx')
    
    filehead = {'Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('ERA Toolbox v%s',REL.eraver);
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',REL.filename);
    filehead{end+1} = sprintf('Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(stddevtable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',REL.eraver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Chains: %d, Iterations: %d',...
        REL.nchains,REL.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s', strcat('Label,Beteen-Person Std Dev',...
        ',Within-Person Std Dev,ICC'));
    fprintf(fid,' \n');
    
    for i = 1:height(stddevtable)
         formatspec = '%s,%s,%s,%s\n';
         fprintf(fid,formatspec,char(stddevtable{i,1}),...
             char(stddevtable{i,2}), char(stddevtable{i,3}),...
             char(stddevtable{i,4}));
    end
    
    fclose(fid);
    
end


end

% 
% function depout = era_dep(var_u, var_e, obs)
% %calculate dependability from a Stanfit object 
% 
% depout = var_u.^2 ./ (var_u.^2 + (var_e.^2./obs));
% 
% end
% 
function iccout = era_icc(var_u,var_e)
%calculate iccs from a Stanfit object

iccout = var_u.^2 ./ (var_u.^2 + var_e.^2);

end
