function [era_data,relerr] = era_relsummary(varargin)
%
%Summarize reliabilty of era_data
%
%era_relsummary('era_data',era_data,'depcutoff',depcutoff,...
%   'meascutoff',meascutoff,'depcentmeas',depcentmeas)
%
%Last Modified 8/18/16
%
%Inputs
% era_data - ERA Toolbox data structure array. Variance components should
%  be included.
% depcutoff - dependabiltiy threshold for considering data reliable
% meascutoff - which estimate of dependability to use for cutoff (1 - lower
%  limit of confidence interval, 2 - point estimate, 3 - upper limit of
%  confidence interval)
% depcentmeas - which measure of central tendency to use for estimating
%  overall score dependability after applying trial cutoffs
%
%Optional Input
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% era_data - rel field will be added to era_data that includes a summary of
%  reliability information for the data
% relerr - structure array with information about reliabitliy analyses.
%  Warnings will be provided to the user if they apply to the analyzed
%  dataset
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
% by Peter Clayson (7/24/16)
% peter.clayson@gmail.com
%
%7/30/16 PC
% bug fix: crashing when running analysis 1. Data were not being properly
%  indexed.
%
%8/14/16 PC
% bug fix: data table not being properly indexed when only analyzing one
%  event 


%somersault through inputs
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_dep for more information about inputs'));
    end
    
    %check if era_data was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('era_data',varargin),1);
    if ~isempty(ind)
        era_data = varargin{ind+1}; 
    else 
        error('varargin:era_data',... %Error code and associated error
            strcat('WARNING: era_data not specified \n\n',... 
            'Please input era_data (ERA Toolbox data structure array).\n',...
            'See help era_relsummary for more information \n'));
    end
    
    %check if depcutoff was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('depcutoff',varargin),1);
    if ~isempty(ind)
        depcutoff = varargin{ind+1}; 
    else 
        error('varargin:depvalue',... %Error code and associated error
            strcat('WARNING: depvalue not specified \n\n',... 
            'Please input the dependability threshold for trial cutoffs.\n',...
            'See help era_relsummary for more information \n'));
    end
    
    %check if meascutoff was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('meascutoff',varargin),1);
    if ~isempty(ind)
        meascutoff = varargin{ind+1}; 
    else 
         error('varargin:meascutoff',... %Error code and associated error
            strcat('WARNING: Which dependability estimate not specified \n\n',... 
            'Please enter a valid input for meascutoff\n',...
            '1 - lower limit of credible interval\n',...
            '2 - point estimate of credible interval\n',...
            '3 - upper limit of credible interval\n',...
            'See help era_relsummary for more information \n'));
    end
    
    %check if depcentmeas was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('depcentmeas',varargin),1);
    if ~isempty(ind)
        depcentmeas = varargin{ind+1}; 
    else 
         error('varargin:depcentmeas',... %Error code and associated error
            strcat('WARNING: Which central tendency to use for ',... 
            'estimating overall score dependability\n',...
            'Please ented a valid input for depcentmeas\n',...
            '1 - mean\n',...
            '2 - median\n',...
            'See help era_relsummary for more information \n'));
    end
    
    %check if CI was specified. 
    %If it is not found, set default as 95%
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help era_depvtrialsplot for more information \n'));
        end
    else 
        ciperc = .95; %default: 95%
    end
end

%place era_data.data in REL to work with
REL = era_data.rel;

%create a data structure for storing outputs
data = struct;

%create variable to specify whether there are bad/unreliable data
%default state: 0, will be changed to 1 if there is a problem
relerr = struct();
relerr.trlcutoff = 0;
relerr.trlmax = 0;

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
        relsummary.meascutoff = strcat('Lower Limit of ',...
            sprintf(' %2.0f',ciperc*100),'% Credible Interval');
    case 2
        relsummary.meascutoff = 'Point Estimate';
    case 3
        relsummary.meascutoff = strcat('Upper Limit of ',...
            sprintf(' %2.0f',ciperc*100),'% Credible Interval');
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
                    'obs',[1 ntrials],'CI',ciperc);

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
            
            relerr.trlcutoff = 1;
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
            datatrls = REL.data{1};

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
            datatable = REL.data{1};

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
            
            
            [relsummary.group(gloc).event(eloc).dep.ll,...
                relsummary.group(gloc).event(eloc).dep.m,...
                relsummary.group(gloc).event(eloc).dep.ul] =...
                era_dep('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,...
                'obs',depcent,'CI',ciperc);
            
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

            [relsummary.group(gloc).event(eloc).dep.ll,...
                relsummary.group(gloc).event(eloc).dep.m,...
                relsummary.group(gloc).event(eloc).dep.ul] =...
                era_dep('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,...
                'obs',depcent,'CI',ciperc);
            
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
                
                relerr.trlmax = 1;
                relsummary.group(gloc).event(eloc).rel.m = -1;
                relsummary.group(gloc).event(eloc).rel.ll = -1;
                relsummary.group(gloc).event(eloc).rel.ul = -1;

            end
            
            
        end
        
        %calculate iccs, between-subject standard deviations, and
        %within-subject standard deviations
        [relsummary.group(gloc).event(eloc).icc.ll,...
            relsummary.group(gloc).event(eloc).icc.m,...
            relsummary.group(gloc).event(eloc).icc.ul] = ...
            era_icc('bp',data.g(gloc).e(eloc).sig_u.raw,...
            'wp',data.g(gloc).e(eloc).sig_e.raw,'CI',ciperc);

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
                    'obs',[1 ntrials],'CI',ciperc);

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

                relerr.trlcutoff = 1;
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
                datatrls = REL.data{1};
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
                datatrls = REL.data{1};
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

            datatable = REL.data{1};
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
            
            [relsummary.group(gloc).event(eloc).dep.ll,...
                relsummary.group(gloc).event(eloc).dep.m,...
                relsummary.group(gloc).event(eloc).dep.ul] =...
                era_dep('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,...
                'obs',depcent,'CI',ciperc);
            
            relsummary.group(gloc).event(eloc).dep.meas = depcent;            
            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            
            %just in case a cutoff was extrapolated, the dependability
            %information will be overwritten
            if  relsummary.group(gloc).event(eloc).trlcutoff >...
                    relsummary.group(gloc).event(eloc).trlinfo.max
                
                relerr.trlmax = 1;
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
            
            [relsummary.group(gloc).event(eloc).icc.ll,...
                relsummary.group(gloc).event(eloc).icc.m,...
                relsummary.group(gloc).event(eloc).icc.ul] = ...
                era_icc('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,'CI',ciperc);
            
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
                    'obs',[1 ntrials],'CI',ciperc);

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

                relerr.trlcutoff = 1;
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

            [relsummary.group(gloc).event(eloc).dep.ll,...
                relsummary.group(gloc).event(eloc).dep.m,...
                relsummary.group(gloc).event(eloc).dep.ul] =...
                era_dep('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,...
                'obs',depcent,'CI',ciperc);
            
            relsummary.group(gloc).event(eloc).dep.meas = depcent;

            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
            
            %just in case a cutoff was extrapolated, the dependability
            %information will be overwritten
            if  relsummary.group(gloc).event(eloc).trlcutoff >...
                    relsummary.group(gloc).event(eloc).trlinfo.max

                relerr.trlmax = 1;
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
            [relsummary.group(gloc).event(eloc).icc.ll,...
                relsummary.group(gloc).event(eloc).icc.m,...
                relsummary.group(gloc).event(eloc).icc.ul] = ...
                era_icc('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,'CI',ciperc);

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
                    'obs',[1 ntrials],'CI',ciperc);
                
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
                    relerr.trlcutoff = 1;
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

            [relsummary.group(gloc).event(eloc).dep.ll,...
                relsummary.group(gloc).event(eloc).dep.m,...
                relsummary.group(gloc).event(eloc).dep.ul] =...
                era_dep('bp',data.g(gloc).e(eloc).sig_u.raw,...
                'wp',data.g(gloc).e(eloc).sig_e.raw,...
                'obs',depcent,'CI',ciperc);
            
            relsummary.group(gloc).event(eloc).dep.meas = depcent;                
                relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                
                %check if trial cutoff for an event exceeded the total
                %number of trials for any subjects
                if  relsummary.group(gloc).event(eloc).trlcutoff >...
                        relsummary.group(gloc).event(eloc).trlinfo.max
                    
                    relerr.trlmax = 1;
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
                [relsummary.group(gloc).event(eloc).icc.ll,...
                    relsummary.group(gloc).event(eloc).icc.m,...
                    relsummary.group(gloc).event(eloc).icc.ul] = ...
                    era_icc('bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',data.g(gloc).e(eloc).sig_e.raw,'CI',ciperc);
               
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
era_data.relsummary = relsummary;
era_data.relsummary.data = data;

end

