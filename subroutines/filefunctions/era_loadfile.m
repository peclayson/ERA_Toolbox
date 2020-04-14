function dataout = era_loadfile(varargin)
%Loads and prepares the data file for dependability analyses
%
%era_loadfile('file',filename)
%
%Last Modified 2/25/19
%
%Required Inputs:
% era_prefs - preferences from ERA Toolbox
% era_data - data from ERA Toolbox
%       OR
% file - location of file to be loaded and prepared for dependability
%  analyses
%
%Optional Inputs:
% idcol - column label for the participant id variable (default: 'id')
% meascol - column label for the measurement to be analyzed (default:
%  'measurement')
% groupcol - column label for the group variable. If no label is provided 
%  it is assumed that there is only one group in the data file.
% whichgroups - cell array of labels for the group types found in the 
%  group column to use
% eventcol - column label for the event variable. If no label is provided
%  it is assumed there is only one event type in the data file.
% whichevents - cell array of labels for the event types found in the 
%  event column to use
% timecol - column label for the occasion variable. If no label is provided
%  it is assumed there is only one occasion in the data file.
% whichtimess - cell array of labels for the occasions found in the 
%  occasion column to use
% dataraw - raw data table outputted from era_startproc (so Matlab doesn't
%  have to re-load entire table)
%
%Output:
% dataout - matlab table with prepared data for reliability analysis.
%  Depending on specifications, the table will have 2 to 5 columns.
%  id: Subject ID (string variable)
%  meas: Measurement
%  group: Group (only when specified, string variable)
%  event: Event Type (only when specified, string variable)
%  time: Occasion (only when specified, string variable)

% Copyright (C) 2016-2020 Peter E. Clayson
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
% by Peter Clayson (4/17/16)
% peter.clayson@gmail.com
%
%7/21/16 PC
% Changes associated with new era_prefs and era_data structure arrays
%
%7/30/16 PC
% Added error catch to verify that the data in the measurement column is
%  numeric
%
%1/19/17 PC
% updated copyright
%
%8/13/17 PC
% added error check: verifies that each participant has a measurement for
%  each event type
%
%8/25/17 PC
% changes associated with adding feature to select only a subset of groups
%  and/or events to process
%
%9/25/17 PC
% fixed bug when selecting a subset of groups and/or events based on
%  numerical inputs
%
%10/16/17 PC
% ran into a bug when indexing a whichgroups or whichevents when containing
%  charcter arrays
%
%4/19/18 PC
% see issue #17 on github. fixed bug when indexing whichgroups and 
%  whichevents that was not loaded as cell array 
%
%4/20/18 PC
% see issue #17 on github and comment immediately before this. Making one
%  change to try and fix bug.
%
%10/29/18 PC
% continued issues loading events/groups that don't follow typical formats.
%
%2/25/19 PC
% changes associated with adding functionality for test-retest reliability 

%try to load era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%if era_prefs and era_data were defined, pull inputs
if ~isempty(era_prefs) && ~isempty(era_data)
    file = era_data.raw.filename;
    idcolname = era_data.proc.idheader;
    meascolname = era_data.proc.measheader;
    eventcolname = era_data.proc.eventheader;
    whichevents = era_data.proc.whichevents;
    groupcolname = era_data.proc.groupheader;
    whichgroups = era_data.proc.whichgroups;
    timecolname = era_data.proc.timeheader;
    whichtimes = era_data.proc.whichtimes;
    dataraw = era_data.raw.data; 
end
    

%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if ~isempty(varargin) && (isempty(era_prefs) || isempty(era_data))
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_loadfile for more information on inputs'));
    end
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = cell2mat(varargin(ind+1)); 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('Error: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
   
    %check if id column was specified 
    %If it is not found, use default 'id'
    ind = find(strcmp('idcol',varargin),1);
    if ~isempty(ind)
        idcolname = cell2mat(varargin(ind+1)); 
    else 
        idcolname = 'id';
    end
    
    %check if group column was specified 
    %If it is not found, assume only one group is present in file
    ind = find(strcmp('groupcol',varargin),1);
    if ~isempty(ind)
        groupcolname = cell2mat(varargin(ind+1)); 
    else 
        groupcolname = '';
    end
    
    %check if whichgroups was specified 
    %If it is not found, assume only one group is present in file
    ind = find(strcmp('whichgroups',varargin),1);
    if ~isempty(ind)
        whichgroups = cell2mat(varargin(ind+1)); 
    else 
        whichgroups = '';
    end
    
    %check if measurement column was specified 
    %If it is not found, assume default 'measurement'
    ind = find(strcmp('meascol',varargin),1);
    if ~isempty(ind)
        meascolname = cell2mat(varargin(ind+1)); 
    else 
        meascolname = 'measurement';
    end
    
    %check if event type column was specified 
    %If it is not found, assume only one event type is present in file
    ind = find(strcmp('eventcol',varargin),1);
    if ~isempty(ind)
        eventcolname = cell2mat(varargin(ind+1)); 
    else 
        eventcolname = '';
    end
    
    %check if whichgroups was specified 
    %If it is not found, assume only one group is present in file
    ind = find(strcmp('whichevents',varargin),1);
    if ~isempty(ind)
        whichevents = cell2mat(varargin(ind+1)); 
    else 
        whichevents = '';
    end
    
    %check if occasion type column was specified 
    %If it is not found, assume only one occasion is present in file
    ind = find(strcmp('eventcol',varargin),1);
    if ~isempty(ind)
        timecolname = cell2mat(varargin(ind+1)); 
    else 
        timecolname = '';
    end
    
    %check if whichtimes was specified 
    %If it is not found, assume only one occasion is present in file
    ind = find(strcmp('whichtimes',varargin),1);
    if ~isempty(ind)
        whichtimes = cell2mat(varargin(ind+1)); 
    else 
        whichtimes = '';
    end
    
    %check if dataraw was specified 
    %If it is not found, create an empty variable
    ind = find(strcmp('dataraw',varargin),1);
    if ~isempty(ind)
        dataraw = varargin{ind+1}; 
    else 
        dataraw = '';
    end
    
elseif isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('Error: Optional inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help era_loadfile for more information on inputs'));
    
end %if ~isempty(varargin)

%load file if it has not been done already
if isempty(dataraw)
    dataraw = era_readtable('file',file);
end

%grab the filename to store
[~,filename] = fileparts(file); 

%make sure all of the necessary headers are present in the file then load
%the data into a table to be outputted for analysis
colnames = dataraw.Properties.VariableNames;
dataout = table;
dataout.Properties.Description = filename;

%ensure that the id column is present
if ~sum(strcmpi(colnames,idcolname)) 
    if ~exist('headererror','var')
        headerror{1} = 'Subject ID';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Subject ID';
    end
elseif sum(strcmpi(colnames,idcolname)) 
    dataout.id = dataraw{:,strcmpi(colnames,idcolname)};
    colnames = dataraw.Properties.VariableNames;
end

%ensure that the measurement column is present
if ~sum(strcmpi(colnames,meascolname)) 
    if ~exist('headererror','var')
        headerror{1} = 'Measurement';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Measurement';
    end
elseif sum(strcmpi(colnames,meascolname)) 
    dataout.meas = dataraw{:,strcmpi(colnames,meascolname)};
end

%ensure that the group column is present if the group header should be
%there
if ~sum(strcmpi(colnames,groupcolname)) && ~isempty(groupcolname)
    if ~exist('headererror','var')
        headerror{1} = 'Group';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Group';
    end
elseif sum(strcmpi(colnames,groupcolname)) && ~isempty(groupcolname)
    dataout.group = dataraw{:,strcmpi(colnames,groupcolname)};
end

%ensure that the event column is present if the event header should be
%there
if ~sum(strcmpi(colnames,eventcolname)) && ~isempty(eventcolname)
    if ~exist('headererror','var')
        headerror{1} = 'Event';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Event';
    end
elseif sum(strcmpi(colnames,eventcolname)) && ~isempty(eventcolname)
    dataout.event = dataraw{:,strcmpi(colnames,eventcolname)};
end

%ensure that the time column is present if the time header should be
%there
if ~sum(strcmpi(colnames,timecolname)) && ~isempty(timecolname)
    if ~exist('headererror','var')
        headerror{1} = 'Occasion';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Occasion';
    end
elseif sum(strcmpi(colnames,timecolname)) && ~isempty(timecolname)
    dataout.time = dataraw{:,strcmpi(colnames,timecolname)};
end

%error catch in case headers for the columns needed are not specified. Let
%the user know which columns were problematic
if exist('headerror','var')
    error('varargin:colheaders',... %Error code and associated error
    strcat('Error: Column headers not properly specified \n\n',... 
    'Please specify the headers for\n',char(strjoin(headerror,', ')),'\n',...
    'See help era_loadfile \n'));
end

%verify that the data in the measurement column is numeric
if ~isnumeric(dataout.meas(:))
    error('meas:notnumeric',... %Error code and associated error
    strcat('Error: The data in the measurement column is not numeric\n',...
    'Was the wrong column specified as mesurement?\n\n',...
    'The column specified was',[' ' meascolname],'\n',...
    'Please specify a column with numeric data (ERP measurements)\n',...
    'See help era_loadfile for more information\n'));
end

%if groups, event and/or time columns are specified, check whether only 
%certain groups or event should be processed. If there are specific groups or
%events to process, make sure those groups and events exist in the data (in
%cases which function was not called by gui).
if ~isempty(groupcolname) && ~isempty(era_data.proc.whichgroups)
    for ii = 1:length(era_data.proc.whichgroups)
        if iscell(era_data.proc.whichgroups)
            if ~isnumeric(era_data.proc.whichgroups{ii})
                if ~any(strcmpi(dataout.group,era_data.proc.whichgroups{ii}))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.group==era_data.proc.whichgroups{ii}))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            end
        else
            if ~isnumeric(era_data.proc.whichgroups(ii))
                if ~any(strcmpi(dataout.group,era_data.proc.whichgroups(ii)))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.group==era_data.proc.whichgroups(ii)))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            end
            
        end
    end
    
    try
        if iscell(era_data.proc.whichgroups)
            if ~isnumeric(era_data.proc.whichgroups{:})
                dataout = dataout(ismember(...
                    dataout.group,era_data.proc.whichgroups),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.group,era_data.proc.whichgroups{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.group,[era_data.proc.whichgroups{:}]),:);
                end
            end
        else
            if ~isnumeric(era_data.proc.whichgroups(:))
                dataout = dataout(ismember(...
                    dataout.group,era_data.proc.whichgroups),:);
            else
                dataout = dataout(ismember(...
                    dataout.group,era_data.proc.whichgroups(:)),:);
            end
        end
    catch
        if iscell(era_data.proc.whichgroups)
            if ~isnumeric(era_data.proc.whichgroups{1})
                dataout = dataout(ismember(...
                    dataout.group,era_data.proc.whichgroups),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.group,era_data.proc.whichgroups{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.group,[era_data.proc.whichgroups{:}]),:);
                end
            end
        else
            if ~isnumeric(era_data.proc.whichgroups(1))
                dataout = dataout(ismember(...
                    dataout.group,era_data.proc.whichgroups),:);
            else
                dataout = dataout(ismember(...
                    dataout.group,era_data.proc.whichgroups(:)),:);
            end
        end
        
    end
end

if ~isempty(eventcolname) && ~isempty(era_data.proc.whichevents)
    for ii = 1:length(era_data.proc.whichevents)
        if iscell(era_data.proc.whichevents)
            if ~isnumeric(era_data.proc.whichevents{ii})
                if ~any(strcmpi(dataout.event,era_data.proc.whichevents{ii}))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.event==era_data.proc.whichevents{ii}))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            end
        else
            if ~isnumeric(era_data.proc.whichevents(ii))
                if ~any(strcmpi(dataout.event,era_data.proc.whichevents(ii)))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.event==era_data.proc.whichevents(ii)))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            end
        end
    end
    
    
    
    try
        if iscell(era_data.proc.whichevents)
            if ~isnumeric(era_data.proc.whichevents{:})
                dataout = dataout(ismember(...
                    dataout.event,era_data.proc.whichevents),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.event,era_data.proc.whichevents{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.event,[era_data.proc.whichevents{:}]),:);
                end
            end
        else
            if ~isnumeric(era_data.proc.whichevents(:))
                dataout = dataout(ismember(...
                    dataout.event,era_data.proc.whichevents),:);
            else
                dataout = dataout(ismember(...
                    dataout.event,era_data.proc.whichevents(:)),:);
            end
        end
        
    catch
        if iscell(era_data.proc.whichevents)
            if ~isnumeric(era_data.proc.whichevents{1})
                dataout = dataout(ismember(...
                    dataout.event,era_data.proc.whichevents),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.event,era_data.proc.whichevents{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.event,[era_data.proc.whichevents{:}]),:);
                end
            end
        else
            if ~isnumeric(era_data.proc.whichevents(1))
                dataout = dataout(ismember(...
                    dataout.event,era_data.proc.whichevents),:);
            else
                dataout = dataout(ismember(...
                    dataout.event,era_data.proc.whichevents(:)),:);
            end
        end
    end

end

if ~isempty(timecolname) && ~isempty(era_data.proc.whichtimes)
    for ii = 1:length(era_data.proc.whichtimes)
        if iscell(era_data.proc.whichtimes)
            if ~isnumeric(era_data.proc.whichtimes{ii})
                if ~any(strcmpi(dataout.time,era_data.proc.whichtimes{ii}))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.time==era_data.proc.whichtimes{ii}))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            end
        else
            if ~isnumeric(era_data.proc.whichtimes(ii))
                if ~any(strcmpi(dataout.time,era_data.proc.whichtimes(ii)))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified ocassions to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.time==era_data.proc.whichtimes(ii)))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help era_loadfile for more information\n'));
                end
            end
        end
    end
    
    
    
    try
        if iscell(era_data.proc.whichtimes)
            if ~isnumeric(era_data.proc.whichtimes{:})
                dataout = dataout(ismember(...
                    dataout.time,era_data.proc.whichtimes),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.time,era_data.proc.whichtimes{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.time,[era_data.proc.whichtimes{:}]),:);
                end
            end
        else
            if ~isnumeric(era_data.proc.whichtimes(:))
                dataout = dataout(ismember(...
                    dataout.time,era_data.proc.whichtimes),:);
            else
                dataout = dataout(ismember(...
                    dataout.time,era_data.proc.whichtimes(:)),:);
            end
        end
        
    catch
        if iscell(era_data.proc.whichtimes)
            if ~isnumeric(era_data.proc.whichtimes{1})
                dataout = dataout(ismember(...
                    dataout.time,era_data.proc.whichtimes),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.time,era_data.proc.whichtimes{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.time,[era_data.proc.whichtimes{:}]),:);
                end
            end
        else
            if ~isnumeric(era_data.proc.whichtimes(1))
                dataout = dataout(ismember(...
                    dataout.time,era_data.proc.whichtimes),:);
            else
                dataout = dataout(ismember(...
                    dataout.time,era_data.proc.whichtimes(:)),:);
            end
        end
    end

end

%verify that each participant has at least one measurement per event type
%(this works because it's already been verified that there are no empty
%cells)
if ~isempty(eventcolname)
    datacheck = varfun(@mean,dataout,'InputVariables','meas',...
       'GroupingVariables',{'id','event'});
    eventcount = varfun(@numel,datacheck,'InputVariables','event',...
        'GroupingVariables','id');
    if length(unique(eventcount.GroupCount)) > 1
        eventcount.GroupCount = [];
        eventcount.Properties.VariableNames{2} = 'Number_of_Events';
        display(eventcount);
        error('meas:mismatchedevents',... %Error code and associated error
        strcat('Error: All participants do not have at least one\n',...
        'measurement per event\n',...
        'The ids and number of events found for each participant are printed',...
        '\nabove to help in finding the problem\n',...
        'See the ''Preparing Data'' section of the User Manual for more information\n'));
    end
end

%verify that each participant has at least one measurement per occasion and
%per event
%(this works because it's already been verified that there are no empty
%cells)
if ~isempty(timecolname) && ~isempty(eventcolname)
    datacheck = varfun(@mean,dataout,'InputVariables','meas',...
       'GroupingVariables',{'id','event','time'});
    timecount = varfun(@numel,datacheck,'InputVariables','time',...
        'GroupingVariables','id');
    if length(unique(timecount.GroupCount)) > 1
        timecount.GroupCount = [];
        timecount.Properties.VariableNames{2} = 'Number_of_Events';
        display(timecount);
        error('meas:mismatchedoccasions',... %Error code and associated error
        strcat('Error: All participants do not have at least one\n',...
        'measurement per event per occasion\n',...
        'The ids and number of occasions found for each participant are printed',...
        '\nabove to help in finding the problem\n',...
        'See the ''Preparing Data'' section of the User Manual for more information\n'));
    end
elseif ~isempty(timecolname) && isempty(eventcolname)
        datacheck = varfun(@mean,dataout,'InputVariables','meas',...
       'GroupingVariables',{'id','time'});
    timecount = varfun(@numel,datacheck,'InputVariables','time',...
        'GroupingVariables','id');
    if length(unique(timecount.GroupCount)) > 1
        timecount.GroupCount = [];
        timecount.Properties.VariableNames{2} = 'Number_of_Events';
        display(timecount);
        error('meas:mismatchedoccasions',... %Error code and associated error
        strcat('Error: All participants do not have at least one\n',...
        'measurement per event per occasion\n',...
        'The ids and number of occasions found for each participant are printed',...
        '\nabove to help in finding the problem\n',...
        'See the ''Preparing Data'' section of the User Manual for more information\n'));
    end
end

%turn all of the ids, groups, events, and occasions into strings. ERA  
%toolbox will use this later
if isnumeric(dataout.id(:))
    newid = cellstr(num2str(dataout.id(:)));
    dataout.id = [];
    dataout.id = newid(:);
end

if ~isempty(groupcolname) && isnumeric(dataout.group(:))
    newgroup = cellstr(num2str(dataout.group(:)));
    dataout.group = [];
    dataout.group = newgroup(:);
end

if ~isempty(eventcolname) && isnumeric(dataout.event(:)) 
    newevent = cellstr(num2str(dataout.event(:)));
    dataout.event = [];
    dataout.event = newevent(:);
end

if ~isempty(timecolname) && isnumeric(dataout.time(:)) 
    newtime = cellstr(num2str(dataout.time(:)));
    dataout.time = [];
    dataout.time = newtime(:);
end

end