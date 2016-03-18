function era_relfigures(varargin)
%Creates various figures and tables for dependability data. See user manual
% for more specific information about each figure.
%
%era_relfigures('data',ERAData)
%
%Note: The Statistics and Machine Learning Toolbox is required
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
% by Peter Clayson (3/5/15)
% peter.clayson@gmail.com
%

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
    
    %check if depmeas is provided
    ind = find(strcmp('depmeas',varargin),1);
    if ~isempty(ind)
        depmeas = cell2mat(varargin{ind+1}); 
    else 
        depmeas = 'mean'; %default is 1
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
    ind = find(strcmp('depplotntrials',varargin),1);
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

elseif ~isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Optional inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help era_relfigures for more information on optional inputs'));
    
end %if ~isempty(varargin)

%create a data structure for storing outputs
data = struct;

%create variabel to specify whether there are bad/unreliable data
%default state: 0, will be changed to 1 if there is a problem
poorrel = 0;

%check whether any groups exist
if strcmpi(REL.groups,'none')
    ngroups = 1;
    gnames = cellstr(REL.groups);
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

%check whether any events exist
if strcmpi(REL.events,'none')
    nevents = 1;
    enames = cellstr(REL.events);
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
            lblstr = strsplit(REL.out.labels{i},'_');

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

%create an x-axis for the number of observations
ntrials = plotntrials;
x = 1:ntrials;

%create an empty array for storing information into
mrel = zeros(ntrials,0);
% llrel = zeros(ntrials,0);
% ulrel = zeros(ntrials,0);

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
depplot.Position = [125 630 900 450];
fsize = 16;

%extract the data and create the subplots for depplot
for eloc=1:nevents
    for gloc=1:ngroups
        for trial=1:ntrials
            mrel(trial,gloc) = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial));
%             llrel(trial,gloc) = ...
%                 quantile(reliab(data.g(gloc).e(eloc).sig_u.raw,...
%                 data.g(gloc).e(eloc).sig_e.raw,trial),.025);
%             ulrel(trial,gloc) = ...
%                 quantile(reliab(data.g(gloc).e(eloc).sig_u.raw,...
%                 data.g(gloc).e(eloc).sig_e.raw,trial),.975);
        end
    end
    subplot(yplots,xplots,eloc);   %Need to grab color from the subplot
    h = plot(x,mrel);
%     clines = get(h,'Color');
% 
%     hold on
%     for i = 1:length(clines)
%         plot(x,llrel(:,i),'--','Color',clines{i});
%         plot(x,ulrel(:,i),'--','Color',clines{i});
%     end
%     hold off
    
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

%create empty arrays for storing dependability information
ntrials = length(data.g(gloc).e(eloc).sig_u.raw);
mrel = zeros(0,ntrials);
llrel = zeros(0,ntrials);
ulrel = zeros(0,ntrials);

%store the cutoff in relsummary to pass to other functions more easily
relsummary.depcutoff = depcutoff;

switch analysis
    case 1 %no groups or event types to consider
        
        eloc = 1;
        gloc = 1;
        
        relsummary.group(gloc).name = gnames{gloc};
        relsummary.group(gloc).event(eloc).name = 'measure';
        
        %compute dependability information
        for trial=1:ntrials
            mrel(trial) = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial)); 
            llrel(trial) = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial),.025);
            ulrel(trial) = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial),.975);
        end

        %find the number of trials to reach cutoff based on the
        %lower limit of the confidence interval
        trlcutoff = find(llrel >= depcutoff, 1);
        
        if isempty(trlcutoff)
            
            poorrel = 1;
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
        
        if isempty(trlcutoff)
            
            %Get trial information for those that meet cutoff
            datatrls = REL.data;

            trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

            ind2exclude = trltable.GroupCount(:);

            relsummary.group(gloc).event(eloc).eventgoodids = 'none';
            relsummary.group(gloc).event(eloc).eventbadids =...
                trltable.id(ind2exclude);

            relsummary.group(gloc).goodids = 'none';
            relsummary.group(gloc).badids = ...
                relsummary.group(gloc).event(eloc).eventbadids;

            datatable = REL.data;

            badids = table(relsummary.group(gloc).badids);

            trlcdata = innerjoin(datatable, badids,...
                'LeftKeys', 'id', 'RightKeys', 'Var1',...
                'LeftVariables', {'id' 'meas'});

            trltable = varfun(@length,trlcdata,...
                'GroupingVariables',{'id'});

            trlmean = mean(trltable.GroupCount);
            trlmed = median(trltable.GroupCount);

            %calculate dependability
            switch depmeas
                case 'mean'
                    depcent = trlmean;
                case 'median'
                    depcent = trlmed;
            end

            mdep = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent)); 
            lldep = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent),.025);
            uldep = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent),.975);

            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depmeas;

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

            %calculate dependability
            switch depmeas
                case 'mean'
                    depcent = trlmean;
                case 'median'
                    depcent = trlmed;
            end

            mdep = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent)); 
            lldep = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent),.025);
            uldep = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent),.975);

            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depmeas;

            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            relsummary.group(gloc).event(eloc).goodn = height(goodids);
            
        end
        
        relsummary.group(gloc).event(eloc).icc.m = ...
            mean(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
            data.g(gloc).e(eloc).sig_e.raw)); 
        relsummary.group(gloc).event(eloc).icc.ll = ...
            quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
            data.g(gloc).e(eloc).sig_e.raw),.025); 
        relsummary.group(gloc).event(eloc).icc.ul = ...
            quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
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
        
       	eloc = 1;

        for gloc=1:ngroups

            if eloc == 1
                relsummary.group(gloc).name = gnames{gloc};
            end

            relsummary.group(gloc).event(eloc).name = enames{eloc};

            for trial=1:ntrials
                mrel(trial) = ...
                    mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,trial)); 
                llrel(trial) = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,trial),.025);
                ulrel(trial) = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,trial),.975);
            end

            %find the number of trials to reach cutoff based on the
            %lower limit of the confidence interval
            trlcutoff = find(llrel >= depcutoff, 1);
            
            if isempty(trlcutoff)

                poorrel = 1;
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

            if trlcutoff == -1
                
                datatrls = REL.data;
                ind = strcmp(datatrls.group,gnames{gloc});
                datatrls = datatrls(ind,:);

                trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                ind2exclude = trltable.GroupCount(:);

                relsummary.group(gloc).event(eloc).eventgoodids = 'none';
                relsummary.group(gloc).event(eloc).eventbadids =...
                    trltable.id(ind2exclude);   
                
            else
                
                datatrls = REL.data;
                ind = strcmp(datatrls.group,gnames{gloc});
                datatrls = datatrls(ind,:);

                trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                ind2include = trltable.GroupCount >= trlcutoff;
                ind2exclude = trltable.GroupCount < trlcutoff;

                relsummary.group(gloc).event(eloc).eventgoodids =...
                    trltable.id(ind2include);
                relsummary.group(gloc).event(eloc).eventbadids =...
                    trltable.id(ind2exclude);
                
            end

            
        end

        %find the ids that have enough trials for each event type
        
        for gloc=1:ngroups
        
            relsummary.group(gloc).goodids = ...
                relsummary.group(gloc).event(eloc).eventgoodids;
            relsummary.group(gloc).badids = ...
                relsummary.group(gloc).event(eloc).eventbadids;
            
        end

        for gloc=1:ngroups

            datatable = REL.data;
            ind = strcmp(datatable.group,gnames{gloc});
            datasubset = datatable(ind,:);
            
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

            %calculate dependability
            switch depmeas
                case 'mean'
                    depcent = trlmean;
                case 'median'
                    depcent = trlmed;
            end

            mdep = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent)); 
            lldep = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent),.025);
            uldep = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,depcent),.975);

            relsummary.group(gloc).event(eloc).dep.m = mdep;
            relsummary.group(gloc).event(eloc).dep.ll = lldep;
            relsummary.group(gloc).event(eloc).dep.ul = uldep;
            relsummary.group(gloc).event(eloc).dep.meas = depmeas;
            
            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            
            if ~strcmp(relsummary.group(gloc).goodids,'none')
                relsummary.group(gloc).event(eloc).goodn = height(goodids);
            elseif strcmp(relsummary.group(gloc).goodids,'none')
                relsummary.group(gloc).event(eloc).goodn = 0;
            end
            
            
            relsummary.group(gloc).event(eloc).icc.m = ...
                mean(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw)); 
            relsummary.group(gloc).event(eloc).icc.ll = ...
                quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw),.025); 
            relsummary.group(gloc).event(eloc).icc.ul = ...
                quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
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
        
        gloc = 1;
        for eloc=1:nevents
                
                if eloc == 1
                    relsummary.group(gloc).name = gnames{gloc};
                end
               
                relsummary.group(gloc).event(eloc).name = enames{eloc};
                
                for trial=1:ntrials
                    mrel(trial) = ...
                        mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial)); 
                    llrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.025);
                    ulrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.975);
                end
                
                %find the number of trials to reach cutoff based on the
                %lower limit of the confidence interval
                trlcutoff = find(llrel >= depcutoff, 1);
                
                if isempty(trlcutoff)
            
                    poorrel = 1;
                    trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).rel.m = -1;
                    relsummary.group(gloc).event(eloc).rel.ll = -1;
                    relsummary.group(gloc).event(eloc).rel.ul = -1;
                    
                    datatrls = REL.data{eloc};
                
                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                    ind2include = 'none';
                    ind2exclude = trltable.GroupCount(:);

                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        'none';
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id(ind2exclude);
                    
                elseif ~isempty(trlcutoff)

                    %store information about cutoffs
                    relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                    relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);

                    datatrls = REL.data{eloc};
                
                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                    ind2include = trltable.GroupCount >= trlcutoff;
                    ind2exclude = trltable.GroupCount < trlcutoff;

                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        trltable.id(ind2include);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id(ind2exclude);
                    
                end 
                
        end
        
        tempids = {};
        badids = [];
        for eloc=1:nevents
            if ~strcmp(relsummary.group(gloc).event(eloc).eventgoodids,'none')
            
                tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;

                if eloc > 1
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
        
        if isempty(tempids)
            tempids{1} = 'none';
        end
        
        relsummary.group(gloc).goodids = tempids{1};
        relsummary.group(gloc).badids = unique(badids);
        
        for eloc=1:nevents

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

                %calculate dependability
                switch depmeas
                    case 'mean'
                        depcent = trlmean;
                    case 'median'
                        depcent = trlmed;
                end

                mdep = ...
                    mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,depcent)); 
                lldep = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,depcent),.025);
                uldep = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,depcent),.975);

                relsummary.group(gloc).event(eloc).dep.m = mdep;
                relsummary.group(gloc).event(eloc).dep.ll = lldep;
                relsummary.group(gloc).event(eloc).dep.ul = uldep;
                relsummary.group(gloc).event(eloc).dep.meas = depmeas;
                
                relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                
                if ~strcmp(relsummary.group(gloc).goodids,'none')
                    relsummary.group(gloc).event(eloc).goodn = height(goodids);
                else
                    relsummary.group(gloc).event(eloc).goodn = 0;
                end
                
                relsummary.group(gloc).event(eloc).icc.m = ...
                    mean(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw)); 
                relsummary.group(gloc).event(eloc).icc.ll = ...
                    quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw),.025); 
                relsummary.group(gloc).event(eloc).icc.ul = ...
                    quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
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
        
        for eloc=1:nevents
            
            for gloc=1:ngroups
                
                if eloc == 1
                    relsummary.group(gloc).name = gnames{gloc};
                end
               
                relsummary.group(gloc).event(eloc).name = enames{eloc};
                
                for trial=1:ntrials
                    mrel(trial) = ...
                        mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial)); 
                    llrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.025);
                    ulrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.975);
                end
                
                %find the number of trials to reach cutoff based on the
                %lower limit of the confidence interval
                trlcutoff = find(llrel >= depcutoff, 1);
                
                if isempty(trlcutoff)
            
                    poorrel = 1;
                    trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).rel.m = -1;
                    relsummary.group(gloc).event(eloc).rel.ll = -1;
                    relsummary.group(gloc).event(eloc).rel.ul = -1;
                    
                    datatrls = REL.data{eloc};
                
                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

                    ind2include = 'none';
                    ind2exclude = trltable.GroupCount(:);

                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        'none';
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id(ind2exclude);
                    
                elseif ~isempty(trlcutoff)

                    %store information about cutoffs
                    relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                    relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);

                    datatrls = REL.data{eloc};
                
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

                    if eloc > 1
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
            
            if isempty(tempids)
                tempids{1} = 'none';
            end
            
            relsummary.group(gloc).goodids = tempids{1};
            relsummary.group(gloc).badids = unique(badids);
            
        end
         
        
        for eloc=1:nevents

            for gloc=1:ngroups
                
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

                %calculate dependability
                switch depmeas
                    case 'mean'
                        depcent = trlmean;
                    case 'median'
                        depcent = trlmed;
                end

                mdep = ...
                    mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,depcent)); 
                lldep = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,depcent),.025);
                uldep = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,depcent),.975);

                relsummary.group(gloc).event(eloc).dep.m = mdep;
                relsummary.group(gloc).event(eloc).dep.ll = lldep;
                relsummary.group(gloc).event(eloc).dep.ul = uldep;
                relsummary.group(gloc).event(eloc).dep.meas = depmeas;
                
                relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                
                if ~strcmp(relsummary.group(gloc).event(eloc).eventgoodids,'none')
                    relsummary.group(gloc).event(eloc).goodn = height(goodids);
                else
                    relsummary.group(gloc).event(eloc).goodn = 0;
                end
                
                relsummary.group(gloc).event(eloc).icc.m = ...
                    mean(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw)); 
                relsummary.group(gloc).event(eloc).icc.ll = ...
                    quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw),.025); 
                relsummary.group(gloc).event(eloc).icc.ul = ...
                    quantile(iccfun(data.g(gloc).e(eloc).sig_u.raw,...
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

RELout = REL;
RELout.relsummary = relsummary;

if picc == 1
    iccplot = figure;
    iccplot.Position = [125 65 900 450];
    miccm = zeros(nevents,ngroups);
    lliccm = zeros(nevents,ngroups);
    uliccm = zeros(nevents,ngroups);
    offsetm = zeros(nevents,ngroups);
    for gloc=1:ngroups
       for eloc=1:nevents
           miccm(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.m;
           lliccm(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.m...
               - relsummary.group(gloc).event(eloc).icc.ll;
           uliccm(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.ul...
               - relsummary.group(gloc).event(eloc).icc.m;
           
           if gloc < median(1:ngroups)
               offsetm(eloc,gloc) = eloc - (.4/ngroups);
           elseif gloc > median(1:ngroups)
               offsetm(eloc,gloc) = eloc + (.4/ngroups);
           elseif gloc == median(1:ngroups)
               offsetm(eloc,gloc) = eloc;
           end
           
       end
    end
    
    iccplot = errorbar(offsetm,miccm,lliccm,uliccm,'Marker','.',...
        'MarkerSize',15,'LineWidth',1);
    ylabel('Intraclass Correlation Coefficient','FontSize',fsize);
    
    for i = 1:length(iccplot)
        iccplot(i).LineStyle = 'none';
    end
    
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
    
    pl = legend(iccplot);
    pl.String = gnames;
    
    iccplot(1).Parent.YAxisLocation = 'right';
    camroll(-90);
    
end

if plotbetstddev == 1
    sdplot = figure;
    sdplot.Position = [225 165 900 450];
    msd = zeros(nevents,ngroups);
    llsd = zeros(nevents,ngroups);
    ulsd = zeros(nevents,ngroups);
    offsetm = zeros(nevents,ngroups);
    for gloc=1:ngroups
       for eloc=1:nevents
           msd(eloc,gloc) = relsummary.group(gloc).event(eloc).betsd.m;
           llsd(eloc,gloc) = relsummary.group(gloc).event(eloc).betsd.m...
               - relsummary.group(gloc).event(eloc).betsd.ll;
           ulsd(eloc,gloc) = relsummary.group(gloc).event(eloc).betsd.ul...
               - relsummary.group(gloc).event(eloc).betsd.m;
           
           if gloc < median(1:ngroups)
               offsetm(eloc,gloc) = eloc - (.4/ngroups);
           elseif gloc > median(1:ngroups)
               offsetm(eloc,gloc) = eloc + (.4/ngroups);
           elseif gloc == median(1:ngroups)
               offsetm(eloc,gloc) = eloc;
           end
           
       end
    end
    
    sdplot = errorbar(offsetm,msd,llsd,ulsd,'Marker','.',...
        'MarkerSize',15,'LineWidth',1);
    ylabel('Between-Person Standard Deviation','FontSize',fsize);
    
    for i = 1:length(sdplot)
        sdplot(i).LineStyle = 'none';
    end
    
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
    
    pl = legend(sdplot);
    pl.String = gnames;
    
    sdplot(1).Parent.YAxisLocation = 'right';
    camroll(-90);
    
end

%reminder....
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

%Create table to display both sets of data
label = {};
trlcutoff = {};
dep = {};
overalldep = {};
mintrl = {};
maxtrl = {};
meantrl = {};
goodn = {};
badn = {};
icc = {};
betsd = {};
witsd = {};

%put data together in tables to display

for gloc=1:ngroups
    for eloc=1:nevents
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
        trlcutoff{end+1} = relsummary.group(gloc).event(eloc).trlcutoff;
        dep{end+1} = sprintf('        %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).rel.m,...
            relsummary.group(gloc).event(eloc).rel.ll,...
            relsummary.group(gloc).event(eloc).rel.ul);
        
        overalldep{end+1} = sprintf('        %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).dep.m,...
            relsummary.group(gloc).event(eloc).dep.ll,...
            relsummary.group(gloc).event(eloc).dep.ul);
        
        mintrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.min;
        maxtrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.max;
        meantrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.mean;
        goodn{end+1} = relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = length(relsummary.group(gloc).badids);
        icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).icc.m,...
            relsummary.group(gloc).event(eloc).icc.ll,...
            relsummary.group(gloc).event(eloc).icc.ul);
        
        betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).betsd.m,...
            relsummary.group(gloc).event(eloc).betsd.ll,...
            relsummary.group(gloc).event(eloc).betsd.ul);
        
        witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            relsummary.group(gloc).event(eloc).witsd.m,...
            relsummary.group(gloc).event(eloc).witsd.ll,...
            relsummary.group(gloc).event(eloc).witsd.ul);
        
    end 
end

stddevtable = table(label',goodn',betsd',witsd');

stddevtable.Properties.VariableNames = {'Label','n_Included',...
    'Between_StdDev','Within_StdDev'};

inctrltable = table(label',trlcutoff',dep');

inctrltable.Properties.VariableNames = {'Label','Trial_Cutoff',...
    'Dependability'};

overalltable = table(label',goodn',badn',overalldep',icc',meantrl',...
    mintrl',maxtrl');

overalltable.Properties.VariableNames = {'Label' ...
    'n_Included' 'n_Excluded' ...
    'Dependability' 'ICC' 'Mean_Num_Trials' 'Min_Num_Trials'...
    'Max_Num_Trials'};



%define parameters for figure size
figwidth = 650;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;
name = ['Point and 95% Interval Estimates for the Between-'...
    'and Within-Person Standard Deviations'];
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
    'Between- and Within-Person Standard Deviations',...
    'Position',[0 row figwidth 25]);          

%Start a table
t = uitable('Parent',era_stddev,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(stddevtable));
set(t,'ColumnName',{'Label' 'n Included' 'Between Std Dev'...
    'Within Std Dev'});
set(t,'ColumnWidth',{200 'auto' 140 140});
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
row = figheight - rowspace*2;

era_inctrl= figure('unit','pix','Visible','off',...
  'position',[1150 700 figwidth figheight],...
  'menub','no',...
  'name',...
  'Results of Increasing Trials the Number of Trials on Dependability',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(era_inctrl,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    sprintf('Dependability Analyses, %0.2f Cutoff',...
    RELout.relsummary.depcutoff),...
    'Position',[0 row figwidth 25]);          

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
    'Dependability' 'ICC' 'Mean # of Trials' 'Min # of Trials'...
    'Max # of Trials'});
set(t,'ColumnWidth',{'auto' 'auto' 'auto' 'auto' 100 'auto' 'auto' 'auto'});
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

if showstddevt == 1
    set(era_stddev,'Visible','on');
end

if pdep == 1
    set(depplot,'Visible','on');
end

if showinct == 1
    set(era_inctrl,'Visible','on');
end

if showoverallt == 1
    set(era_overall,'Visible','on');
end

if poorrel == 1
    
end

end

function era_saveoveralltable(varargin)

REL = varargin{3};
overalltable = varargin{4};

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

if strcmp(ext,'.xlsx')
    
    filehead = {'Dependability Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('Dataset: %s',REL.filename);
    filehead{end+1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(overalltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Dependability Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'\n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'%s', strcat('Label,N Included,N Excluded,',...
        'Dependability,ICC,Mean Num of Trials,Min Num of Trials,',...
        'Max Num of Trials'));
    fprintf(fid,'\n');
    for i = 1:height(overalltable)
         formatspec = '%s,%d,%d,%s,%s,%0.4f,%d,%d\n';
         fprintf(fid,formatspec,char(overalltable{i,1}),...
             cell2mat(overalltable{i,2}),cell2mat(overalltable{i,3}),...
             cell2mat(overalltable{i,4}),char(overalltable{i,5}),...
             cell2mat(overalltable{i,6}),cell2mat(overalltable{i,7}),...
             cell2mat(overalltable{i,8}));
    end
    
    fclose(fid);
    
end

end

function era_saveids(varargin)

REL = varargin{3};

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

if strcmp(ext,'.xlsx')
    
    datap{1,1} = 'Data Generated on';
    datap{2,1} = datestr(clock); 
    datap{3,1} = sprintf('Dataset: %s',REL.filename);
    datap{4,1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    datap{5,1}='';
    datap{6,1}='';
    datap{7,1} = 'Good IDs'; datap{7,2} = 'Bad IDs';
    gids = REL.relsummary.group.goodids;
    for i = 1:length(gids)
        datap{i+7,1}=char(gids(i));
    end
    bids = REL.relsummary.group.badids;
    for i = 1:length(bids)
        datap{i+7,2}=char(bids(i));
    end
    
    xlswrite(fullfile(savepath,savename),datap);
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Data Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'\n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','Good IDs,Bad IDs');
    
    gids = REL.relsummary.group.goodids;
    bids = REL.relsummary.group.badids;

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

REL = varargin{3};
incltrltable = varargin{4};

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

if strcmp(ext,'.xlsx')
    
    filehead = {'Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('Dataset: %s',REL.filename);
    filehead{end+1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(incltrltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'\n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'%s', strcat('Label,Trial Cutoff,Dependability'));
    fprintf(fid,'\n');
    
    for i = 1:height(incltrltable)
         formatspec = '%s,%d,%s\n';
         fprintf(fid,formatspec,char(incltrltable{i,1}),...
             cell2mat(incltrltable{i,2}),char(incltrltable{i,3}));
    end
    
    fclose(fid);
    
end


end


function era_savestddevtable(varargin)

REL = varargin{3};
stddevtable = varargin{4};

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

if strcmp(ext,'.xlsx')
    
    filehead = {'Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('Dataset: %s',REL.filename);
    filehead{end+1} = sprintf('Dependability Cutoff: %0.2f',...
        REL.relsummary.depcutoff);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(stddevtable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'\n');
    fprintf(fid,'Dataset: %s\n',REL.filename);
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        REL.relsummary.depcutoff);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'%s', strcat('Label,n Included,Beteen-Person Std Dev',...
        ',Within-Person Std Dev'));
    fprintf(fid,'\n');
    
    for i = 1:height(stddevtable)
         formatspec = '%s,%d,%s,%s\n';
         fprintf(fid,formatspec,char(stddevtable{i,1}),...
             cell2mat(stddevtable{i,2}),char(stddevtable{i,3}),...
             char(stddevtable{i,4}));
    end
    
    fclose(fid);
    
end


end


function depout = reliab(var_u, var_e, obs)

depout = var_u.^2 ./ (var_u.^2 + (var_e.^2./obs));

end

function iccout = iccfun(var_u,var_e)

iccout = var_u.^2 ./ (var_u.^2 + var_e.^2);

end

function dep = depall(lme,num)

[obsvar,resvar] = covarianceParameters(lme);

dep = cell2mat(obsvar)/(cell2mat(obsvar)+(resvar/num));

end