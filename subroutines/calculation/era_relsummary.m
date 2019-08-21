function [era_data,relerr] = era_relsummary(varargin)
%Summarize reliabilty of era_data
%
%era_relsummary('era_data',era_data,'analysis','sing',...
%   'depcutoff',depcutoff,'meascutoff',meascutoff,...
%   'depcentmeas',depcentmeas)
%
%Last Modified 8/3/19
%
%Inputs
% era_data - ERA Toolbox data structure array. Variance components should
%  be included
% analysis - 'sing' single occasion data; 'trt' multiple occasion data
%
% EITHER (for analysis type 'sing')
%  depcutoff - dependabiltiy threshold for considering data reliable
%  meascutoff - which estimate of dependability to use for cutoff (1 - lower
%   limit of confidence interval, 2 - point estimate, 3 - upper limit of
%   confidence interval)
%  depcentmeas - which measure of central tendency to use for estimating
%   overall score dependability after applying trial cutoffs
%
% OR (for analysis type 'trt')
%  gcoeff - g-theory coefficient to calculate 1 - dependability, 2 -
%   generalizability
%  reltype - reliability coefficient to calculate 1 - equivalence, 2 -
%   stability
%  relcutoff - reliability threshold for considering data reliable
%  meascutoff - which estimate of dependability to use for cutoff (1 - lower
%   limit of confidence interval, 2 - point estimate, 3 - upper limit of
%   confidence interval)
%  relcentmeas - which measure of central tendency to use for estimating
%   overall score reliability after applying trial cutoffs
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

% Copyright (C) 2016-2019 Peter E. Clayson
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
%
%8/21/16 PC
% increased trial estimation for cutoffs from max number of trials +100 to
%  max+1000
%
%10/23/16 PC
% error catch for when no goodids are left even after extrapolating +1000
%  trials out
%
%11/10/16 pc
% bug fix: depending on version of Matlab there was a problem indexing a
%  table for the group analysis when there only two groups
%
%1/19/17 PC
% updated copyright
%
%8/25/17 PC
% minor changes after adding feature to select subsets of groups and events
%  to process
%
%6/22/18 PC
% added standard deviation to trial summaries
%
%6/21/19 PC
% started working on adding trt capabilities
%
%8/2/19 PC
% finished adding trt capability
%
%8/3/19 PC
% now store confidence interval used for reliability estimates
%
%8/21/19 PC
% bug fix: was overwriting overall reliability with trial cutoff
%  reliability for trt data

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
    
    %check if the dependability cutoff was specified
    ind = find(strcmp('analysis',varargin),1);
    if ~isempty(ind)
        relanalysis = varargin{ind+1};
    else
        error('varargin:noanalysistype',... %Error code and associated error
            strcat('WARNING: analysis not specified \n\n',...
            'Please input the analysis to be run: ''sing'' or ''trt'' \n'));
    end
end

switch relanalysis
    case 'sing'
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
    case 'trt'
        %check if gcoeff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('gcoeff',varargin),1);
        if ~isempty(ind)
            gcoeff = varargin{ind+1};
        else
            error('varargin:gcoeff',... %Error code and associated error
                strcat('WARNING: gcoeff not specified \n\n',...
                'Please input the g-theory coefficient to calculate.\n',...
                'See help era_relsummary for more information \n'));
        end
        
        %check if reltype was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('reltype',varargin),1);
        if ~isempty(ind)
            reltype = varargin{ind+1};
        else
            error('varargin:reltype',... %Error code and associated error
                strcat('WARNING: reltype not specified \n\n',...
                'Please input the reliability coefficient to calculate.\n',...
                'See help era_relsummary for more information \n'));
        end
        
        %check if relcutoff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('relcutoff',varargin),1);
        if ~isempty(ind)
            relcutoff = varargin{ind+1};
        else
            error('varargin:relcutoff',... %Error code and associated error
                strcat('WARNING: relcutoff not specified \n\n',...
                'Please input the reliability threshold for trial cutoffs.\n',...
                'See help era_relsummary for more information \n'));
        end
        
        %check if meascutoff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('meascutoff',varargin),1);
        if ~isempty(ind)
            meascutoff = varargin{ind+1};
        else
            error('varargin:meascutoff',... %Error code and associated error
                strcat('WARNING: Which estimate not specified \n\n',...
                'Please enter a valid input for meascutoff\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help era_relsummary for more information \n'));
        end
        
        %check if relcentmeas was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('relcentmeas',varargin),1);
        if ~isempty(ind)
            relcentmeas = varargin{ind+1};
        else
            error('varargin:relcentmeas',... %Error code and associated error
                strcat('WARNING: Which central tendency to use for ',...
                'estimating overall score reliability\n',...
                'Please ented a valid input for relcentmeas\n',...
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


switch relanalysis
    case 'sing'
        %place era_data.data in REL to work with
        REL = era_data.rel;
        
        %create a data structure for storing outputs
        data = struct;
        
        %create variable to specify whether there are bad/unreliable data
        %default state: 0, will be changed to 1 if there is a problem
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;
        
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
        relsummary.ciperc = ciperc;
        
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
                ntrials = max(trltable.GroupCount(:)) + 1000;
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
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(badids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        errordlg(dlg, 'Data do not meet reliability threshold');
                        relerr.nogooddata = 1;
                        return;
                    end
                    
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
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = 0;
                    
                elseif ~isempty(trlcutoff)
                    
                    %Get trial information for those that meet cutoff
                    datatrls = REL.data;
                    
                    try
                        trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                    catch
                        trltable = varfun(@length,datatrls{:},'GroupingVariables',{'id'});
                    end
                    
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
                    
                    if iscell(datatable)
                        datatable = REL.data{1};
                    end
                    
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
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
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
                    ntrials = max(trltable.GroupCount(:)) + 1000;
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
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
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
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        
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
                    
                    try
                        datatable = REL.data{1};
                    catch
                        datatable = REL.data;
                    end
                    ind = strcmp(datatable.group,gnames{gloc});
                    datasubset = datatable(ind,:);
                    
                    %only factor in the trial counts from those with good data
                    
                    if ~strcmp(relsummary.group(gloc).goodids,'none')
                        goodids = table(relsummary.group(gloc).goodids);
                    elseif strcmp(relsummary.group(gloc).goodids,'none')
                        goodids = table(relsummary.group(gloc).badids);
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        errordlg(dlg, 'Data do not meet reliability threshold');
                        relerr.nogooddata = 1;
                        return;
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
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
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
                    ntrials = max(trltable.GroupCount(:)) + 1000;
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
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        errordlg(dlg, 'Data do not meet reliability threshold');
                        relerr.nogooddata = 1;
                        return;
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
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
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
                        ntrials = max(trltable.GroupCount(:)) + 1000;
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
                        
                        %check whether there are not enough good data after applying
                        %reliability threshold and extrapolating
                        if isempty(goodids)
                            dlg = {'Data do not reach reliability threshold after extrapolation';...
                                'Set a lower reliability threshold'};
                            errordlg(dlg, 'Data do not meet reliability threshold');
                            relerr.nogooddata = 1;
                            return;
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
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
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
        
    case 'trt'
        %place era_data.data in REL to work with
        REL = era_data.rel;
        
        %create a data structure for storing outputs
        data = struct;
        
        %create variable to specify whether there are bad/unreliable data
        %default state: 0, will be changed to 1 if there is a problem
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;
        
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
                
                data.g(gloc).e(eloc).label = REL.out.labels{gloc};
                data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
                data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
                data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
                data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
                data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
                data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
                data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
                data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
                data.g(gloc).e(eloc).elabel = cellstr('none');
                data.g(gloc).glabel = gnames(gloc);
                
            case 2 %2 - possible multiple groups but no event types to consider
                
                eloc = 1;
                
                for gloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{gloc};
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
                    data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
                    data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
                    data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
                    data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
                    data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
                    data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
                    data.g(gloc).e(eloc).elabel = cellstr('none');
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 3 %3 - possible event types but no groups to consider
                
                gloc = 1;
                
                for eloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{eloc};
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
                    data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,eloc);
                    data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,eloc);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,eloc);
                    data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,eloc);
                    data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,eloc);
                    data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,eloc);
                    data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,eloc);
                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 4 %4 - possible groups and event types to consider
                for ii=1:length(REL.out.labels)
                    
                    %use the underscores that were added in era_computerel to
                    %differentiate where the group and event the data are for
                    lblstr = strsplit(REL.out.labels{ii},'_;_');
                    
                    eloc = find(ismember(enames,lblstr(2)));
                    gloc = find(ismember(gnames,lblstr(1)));
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{ii};
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,ii);
                    data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,ii);
                    data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,ii);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,ii);
                    data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,ii);
                    data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,ii);
                    data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,ii);
                    data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,ii);
                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
        end %switch analysis
        
        %store reliability info in relsummary to pass to other functions more easily
        relsummary.relcutoff = relcutoff;
        relsummary.reltype = reltype;
        relsummary.ciperc = ciperc;
        
        if reltype == 1
            relsummary.reltype_name = 'ic';
        elseif reltype == 2
            relsummary.reltype_name = 'trt';
        end
        
        relsummary.gcoeff = gcoeff;
        
        if gcoeff == 1
            relsummary.gcoeff_name = 'dep';
        elseif gcoeff == 2
            relsummary.gcoeff_name = 'gen';
        end
        
        
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
        
        
        %compute reliability data for each group and event
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
                    %create empty arrays for storing reliability information
                    trltable = varfun(@length,REL.data{1},'GroupingVariables',{'id'});
                catch
                    %create empty arrays for storing reliability information
                    trltable = varfun(@length,REL.data,'GroupingVariables',{'id'});
                end
                
                %compute reliabiltiy
                ntrials = max(trltable.GroupCount(:)) + 1000;
                [llrel,mrel,ulrel] = era_rel_trt(...
                    'gcoeff',gcoeff,...
                    'reltype',reltype,...
                    'bp',data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[1 ntrials],'CI',ciperc);
                
                %find the number of trials to reach cutoff
                switch meascutoff
                    case 1
                        trlcutoff = find(llrel >= relcutoff, 1);
                    case 2
                        trlcutoff = find(mrel >= relcutoff, 1);
                    case 3
                        trlcutoff = find(ulrel >= relcutoff, 1);
                end
                
                %see whether the trial cutoff was found. If not store all values as
                %-1. If it is found, store the reliability information about the
                %cutoffs.
                if isempty(trlcutoff)
                    
                    relerr.trlcutoff = 1;
                    trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).relcut.m = -1;
                    relsummary.group(gloc).event(eloc).relcut.ll = -1;
                    relsummary.group(gloc).event(eloc).relcut.ul = -1;
                    
                elseif ~isempty(trlcutoff)
                    
                    %store information about cutoffs
                    relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                    relsummary.group(gloc).event(eloc).relcut.m = mrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).relcut.ll = llrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).relcut.ul = ulrel(trlcutoff);
                    
                end
                
                %find the participants without enough trials based on the cutoffs
                if isempty(trlcutoff) %if the cutoff was not found in the data
                    
                    try
                        %Get information for the participants
                        datatrls = REL.data;
                    catch
                        datatrls = REL.data{1};
                    end
                    
                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                    
                    %define all of the data as bad, because the cutoff wasn't even
                    %found
                    relsummary.group(gloc).event(eloc).eventgoodids = 'none';
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id;
                    
                    relsummary.group(gloc).goodids = 'none';
                    relsummary.group(gloc).badids = ...
                        relsummary.group(gloc).event(eloc).eventbadids;
                    
                    %calculate the reliability estimates for the overall data
                    try
                        datatable = REL.data;
                    catch
                        datatable = REL.data{1};
                    end
                    
                    badids = table(relsummary.group(gloc).badids);
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(badids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        errordlg(dlg, 'Data do not meet reliability threshold');
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datatable, badids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        era_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = 0;
                    
                elseif ~isempty(trlcutoff)
                    
                    %Get trial information for those that meet cutoff
                    datatrls = REL.data;
                    
                    try
                        trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                    catch
                        trltable = varfun(@length,datatrls{:},'GroupingVariables',{'id'});
                    end
                    
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
                    
                    if iscell(datatable)
                        datatable = REL.data{1};
                    end
                    
                    goodids = table(relsummary.group(gloc).goodids);
                    
                    trlcdata = innerjoin(datatable, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        era_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    
                    %just in case a cutoff was extrapolated, the reliability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).relcut.m = -1;
                        relsummary.group(gloc).event(eloc).relcut.ll = -1;
                        relsummary.group(gloc).event(eloc).relcut.ul = -1;
                        
                    end
                    
                    %calculate iccs, between-subject standard deviations, and
                    %within-subject standard deviations
                    [relsummary.group(gloc).event(eloc).icc.ll,...
                        relsummary.group(gloc).event(eloc).icc.m,...
                        relsummary.group(gloc).event(eloc).icc.ul] = ...
                        era_icc('bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'wp',data.g(gloc).e(eloc).sig_err.raw,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_id.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,.025);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,.975);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_err.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,.025);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,.975);
                    
                end
                
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
                        %create empty arrays for storing reliability information
                        trltable = varfun(@length,REL.data{1},...
                            'GroupingVariables',{'id'});
                    catch
                        %create empty arrays for storing reliability information
                        trltable = varfun(@length,REL.data,...
                            'GroupingVariables',{'id'});
                    end
                    
                    
                    %compute reliabiltiy
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    [llrel,mrel,ulrel] = era_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',[1 ntrials],'CI',ciperc);
                    
                    %find the number of trials to reach cutoff
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= relcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= relcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= relcutoff, 1);
                    end
                    
                    %see whether the trial cutoff was found. If not store all
                    %values as -1. If it is found, store the reliability
                    %information about the cutoffs.
                    
                    if isempty(trlcutoff)
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).relcut.m = -1;
                        relsummary.group(gloc).event(eloc).relcut.ll = -1;
                        relsummary.group(gloc).event(eloc).relcut.ul = -1;
                        
                    elseif ~isempty(trlcutoff)
                        
                        %store information about cutoffs
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).relcut.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcut.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcut.ul = ulrel(trlcutoff);
                        
                    end
                    
                    %if the trial cutoff was not found in the specified data
                    if trlcutoff == -1
                        
                        %store all the ids as bad
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        
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
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        
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
                
                %calculate the reliability for the overall data for each group
                for gloc=1:ngroups
                    
                    try
                        datatable = REL.data{1};
                    catch
                        datatable = REL.data;
                    end
                    
                    ind = strcmp(datatable.group,gnames{gloc});
                    datasubset = datatable(ind,:);
                    
                    %only factor in the trial counts from those with good data
                    
                    if ~strcmp(relsummary.group(gloc).goodids,'none')
                        goodids = table(relsummary.group(gloc).goodids);
                    elseif strcmp(relsummary.group(gloc).goodids,'none')
                        goodids = table(relsummary.group(gloc).badids);
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        errordlg(dlg, 'Data do not meet reliability threshold');
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datasubset, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        era_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the reliability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).relcut.m = -1;
                        relsummary.group(gloc).event(eloc).relcut.ll = -1;
                        relsummary.group(gloc).event(eloc).relcut.ul = -1;
                        
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
                        era_icc('bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'wp',data.g(gloc).e(eloc).sig_err.raw,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_id.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,.025);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,.975);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_err.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,.025);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,.975);
                    
                    
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
                    
                    
                    %create empty arrays for storing reliability information
                    trltable = varfun(@length,...
                        REL.data(strcmp(REL.data.event,enames{eloc}),:),...
                        'GroupingVariables',{'id'});
                    
                    
                    %compute reliabiltiy
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    [llrel,mrel,ulrel] = era_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',[1 ntrials],'CI',ciperc);
                    
                    %find the number of trials to reach cutoff
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= relcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= relcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= relcutoff, 1);
                    end
                    
                    
                    if isempty(trlcutoff) %if a cutoff wasn't found
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).relcut.m = -1;
                        relsummary.group(gloc).event(eloc).relcut.ll = -1;
                        relsummary.group(gloc).event(eloc).relcut.ul = -1;
                        
                        datatrls = REL.data(strcmp(REL.data.event,enames{eloc}),:);                        
                        
                        %store all of the participant ids as bad, becuase none
                        %reached the cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            'none';
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id;
                        
                    elseif ~isempty(trlcutoff) %if a cutoff was found
                        
                        %store information about cutoffs
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).relcut.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcut.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcut.ul = ulrel(trlcutoff);
                        
                        datatrls = REL.data(strcmp(REL.data.event,enames{eloc}),:);
                        
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
                    
                    %get trial information to use for reliability estimates. only
                    %participants with enough data will contribute toward trial
                    %counts
                    datatable = REL.data(strcmp(REL.data.event,enames{eloc}),:);
                    
                    if ~strcmp(relsummary.group(gloc).goodids,'none')
                        goodids = table(relsummary.group(gloc).goodids);
                    else
                        goodids = table(relsummary.group(gloc).badids);
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        errordlg(dlg, 'Data do not meet reliability threshold');
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datatable, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        era_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the reliability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                        
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
                        era_icc('bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'wp',data.g(gloc).e(eloc).sig_err.raw,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_id.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,.025);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,.975);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_err.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,.025);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,.975);
                    
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
                        
                        %create empty arrays for storing reliability information
                        trltable = varfun(@length,REL.data(all(...
                            strcmp(REL.data.group,gnames{gloc}) &... 
                            strcmp(REL.data.event,enames{eloc})...
                            ,2),:),...
                            'GroupingVariables',{'id'});
                        
                        %compute reliability
                        ntrials = max(trltable.GroupCount(:)) + 1000;
                        [llrel,mrel,ulrel] = era_rel_trt(...
                            'gcoeff',gcoeff,...
                            'reltype',reltype,...
                            'bp',data.g(gloc).e(eloc).sig_id.raw,...
                            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                            'err',data.g(gloc).e(eloc).sig_err.raw,...
                            'obs',[1 ntrials],'CI',ciperc);
                        
                        %find the number of trials to reach cutoff
                        switch meascutoff
                            case 1
                                trlcutoff = find(llrel >= relcutoff, 1);
                            case 2
                                trlcutoff = find(mrel >= relcutoff, 1);
                            case 3
                                trlcutoff = find(ulrel >= relcutoff, 1);
                        end
                        
                        %if a cutoff was not found
                        if isempty(trlcutoff)
                            
                            %store all the participant ids as having bad data
                            relerr.trlcutoff = 1;
                            trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                            
                            datatrls = REL.data(all(...
                                strcmp(REL.data.group,gnames{gloc}) &...
                                strcmp(REL.data.event,enames{eloc})...
                                ,2),:);
                            
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
                            relsummary.group(gloc).event(eloc).relcutoff.m = mrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).relcutoff.ll = llrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).relcutoff.ul = ulrel(trlcutoff);
                            
                            datatrls = REL.data(all(...
                                strcmp(REL.data.group,gnames{gloc}) &...
                                strcmp(REL.data.event,enames{eloc})...
                                ,2),:);
                            
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
                
                %calculate reliability estimates for the overall data for each
                %group and event type
                for eloc=1:nevents
                    
                    for gloc=1:ngroups
                        
                        %pull the data to calculate the trial information from
                        %participants with good data
                        datatable = REL.data(all(...
                            strcmp(REL.data.group,gnames{gloc}) &...
                            strcmp(REL.data.event,enames{eloc})...
                            ,2),:);
                        
                        ind = strcmp(datatable.group,gnames{gloc});
                        datasubset = datatable(ind,:);
                        
                        if ~strcmp(relsummary.group(gloc).goodids,'none')
                            goodids = table(relsummary.group(gloc).goodids);
                        else
                            goodids = table(relsummary.group(gloc).badids);
                        end
                        
                        %check whether there are not enough good data after applying
                        %reliability threshold and extrapolating
                        if isempty(goodids)
                            dlg = {'Data do not reach reliability threshold after extrapolation';...
                                'Set a lower reliability threshold'};
                            errordlg(dlg, 'Data do not meet reliability threshold');
                            relerr.nogooddata = 1;
                            return;
                        end
                        
                        trlcdata = innerjoin(datasubset, goodids,...
                            'LeftKeys', 'id', 'RightKeys', 'Var1',...
                            'LeftVariables', {'id' 'meas'});
                        
                        trltable = varfun(@length,trlcdata,...
                            'GroupingVariables',{'id'});
                        
                        trlmean = mean(trltable.GroupCount);
                        trlmed = median(trltable.GroupCount);
                        
                        %calculate reliability using either the mean or median
                        %(based on user input)
                        switch relcentmeas
                            case 1
                                relcent = trlmean;
                            case 2
                                relcent = trlmed;
                        end
                        
                        [relsummary.group(gloc).event(eloc).rel.ll,...
                            relsummary.group(gloc).event(eloc).rel.m,...
                            relsummary.group(gloc).event(eloc).rel.ul] = ...
                            era_rel_trt(...
                            'gcoeff',gcoeff,...
                            'reltype',reltype,...
                            'bp',data.g(gloc).e(eloc).sig_id.raw,...
                            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                            'err',data.g(gloc).e(eloc).sig_err.raw,...
                            'obs',relcent,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).rel.meas = relcent;
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                        
                        %check if trial cutoff for an event exceeded the total
                        %number of trials for any subjects
                        if  relsummary.group(gloc).event(eloc).trlcutoff >...
                                relsummary.group(gloc).event(eloc).trlinfo.max
                            
                            relerr.trlmax = 1;
                            relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                            
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
                            era_icc('bp',data.g(gloc).e(eloc).sig_id.raw,...
                            'wp',data.g(gloc).e(eloc).sig_err.raw,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).betsd.m = ...
                            mean(data.g(gloc).e(eloc).sig_id.raw);
                        relsummary.group(gloc).event(eloc).betsd.ll = ...
                            quantile(data.g(gloc).e(eloc).sig_id.raw,.025);
                        relsummary.group(gloc).event(eloc).betsd.ul = ...
                            quantile(data.g(gloc).e(eloc).sig_id.raw,.975);
                        
                        relsummary.group(gloc).event(eloc).witsd.m = ...
                            mean(data.g(gloc).e(eloc).sig_err.raw);
                        relsummary.group(gloc).event(eloc).witsd.ll = ...
                            quantile(data.g(gloc).e(eloc).sig_err.raw,.025);
                        relsummary.group(gloc).event(eloc).witsd.ul = ...
                            quantile(data.g(gloc).e(eloc).sig_err.raw,.975);
                        
                    end
                end
                
                
        end %switch analysis
        
end

%store relsummary information
era_data.relsummary = relsummary;
era_data.relsummary.data = data;

end

