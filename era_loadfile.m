function dataout = era_loadfile(varargin)
%Loads and prepares the data file for dependability analyses
%
%era_loadfile('file',filename)
%
%Last Modified 4/18/16
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
% eventcol - column label for the event variable. If no label is provided
%  it is assumed there is only one event type in the data file.
% dataraw - raw data table outputted from era_startproc (so Matlab doesn't
%  have to re-load entire table)
%
%Output:
% dataout - matlab table with prepared data for reliability analysis.
%  Depending on specifications, the table will have 2 to 4 columns.
%  id: Subject ID (string variable)
%  meas: Measurement
%  group: Group (only when specified, string variable)
%  event: Event Type (only when specified, string variable)

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
% by Peter Clayson (4/17/16)
% peter.clayson@gmail.com
%
%7/21/16 PC
% Changes associated with new era_prefs and era_data structure arrays

%try to load era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%if era_prefs and era_data were defined, pull inputs
if ~isempty(era_prefs) && ~isempty(era_data)
    file = era_data.raw.filename;
    idcolname = era_data.proc.idheader;
    groupcolname = era_data.proc.groupheader;
    meascolname = era_data.proc.measheader;
    eventcolname = era_data.proc.eventheader;
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
        strcat('WARNING: File location not specified \n\n',... 
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
    strcat('WARNING: Optional inputs are incomplete \n\n',... 
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

%error catch in case headers for the columns needed are not specified. Let
%the user know which columns were problematic
if exist('headerror','var')
    error('varargin:colheaders',... %Error code and associated error
    strcat('WARNING: Column headers not properly specified \n\n',... 
    'Please specify the headers for\n',char(strjoin(headerror,', ')),'\n',...
    'See help era_loadfile \n'));
end

%turn all of the ids, groups, and events into strings. ERA toolbox will 
%use this later
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

end