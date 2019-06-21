function era_startproc(varargin)
%Initiate Matlab gui to begin processing data in Stan
%
%
%Last Updated 2/25/19
%
%
%Input
% There are no required inputs to execute this script. The user will be
%  asked to identify the location of the data to processed and the columns
%  of interest for analysis
%
%Output
% This script will not output any variables to the workspace. A .mat with
%  processed data from Stan will be saved in a user specified locaiton.
%  This .mat will automatically be loaded into era_startview for viewing 
%  the data
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
% by Peter Clayson (4/18/16)
% peter.clayson@gmail.com
%
%4/20/16 PC
% Add line to print that model converged successfully
% Added ERA Toolbox specific format for saving data (extension: .erat)
%
%7/22/16 PC 
% Pulling out subroutines
%
%7/28/16 PC
% Missed a change related to era_data structure
%
%9/18/16 PC
% Added some tooltips where helpful
%
%10/22/16 PC
% Added an error catch when non-numeric data are specified in the
%  Measurement column
% Fixed error catch for specifying the same column in multiple variables
%
%11/10/16 PC
% bug fix: when trying to rerun on new data when the chains did not
%  converge, ended up with start_proc and start_view both running
%
%1/19/17 PC
% updated copyright
%
%6/24/17 PC
% small change to clarify that a filename needs to be provided for the era
%  output file
%
%8/16/17 PC
% fixed bug when trying to run era_startproc in cmd win without using
%  era_start
% added input to preference window for viewing trace plots prior to saving
%  Stan outputs
%
%8/22/17 PC
% bug fixes for passing input for trace plots to prefs
%
%8/25/17 PC
% new changes to allow user to select a subset of groups and/or events to
%  process
%
%9/25/17 PC
% fixed bug when trying to select a subset of gropus and/or events based on
%  numerical inputs (rather than string)
%
%2/25/19 PC
% adding funcionality to look at test retest reliability
%

%check if era_gui is open. If the user executes era_startproc and skips
%era_start then there will be no gui to close.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    %Grab preferences if they exist
    if ~isempty(varargin)

        %check for era_prefs
        ind = find(strcmp('era_prefs',varargin),1);
        if ~isempty(ind)
            era_prefs = varargin{ind+1}; 
        else 
            %load default preferences for processing and viewing data
            era_prefs = era_defaults;

            %attach the current version number to era_prefs
            era_prefs.ver = era_defineversion;

            %define parameters for figure position
            era_prefs.guis.fsize = get(0,'DefaultTextFontSize');
        end
    end
    
    close(era_gui);
    
elseif isempty(era_gui)
    
    %load default preferences for processing and viewing data
    era_prefs = era_defaults;

    %attach the current version number to era_prefs
    era_prefs.ver = era_defineversion;

    %define parameters for figure position
    era_prefs.guis.fsize = get(0,'DefaultTextFontSize');
    
end

%ask the user to identify the data file to be loaded
[filepart, pathpart] = uigetfile({'*.xlsx;*.xls;*.csv;*.dat;*.txt;*.ods',...
    'All Readable Files'; '*.xlsx','Excel File (.xlsx)';'*.xls',...
    'Excel File 97-2003 (.xls)';'*.csv',...
    'Comma-Separated Value File (.csv)';'*.dat',...
    'Tab-Delimited Text (.dat)';'*.txt',...
    'Raw Text File (.txt)';'*.ods',...
    'OpenDocument Spreadsheet (.ods)'},'Data');

%if the user does not select a file, then take the user back to era_start    
if filepart == 0 
    era_start;
    errordlg('No file selected','File Error');
    return;
end

fprintf('\n\nLoading Data...');
fprintf('\nThis may take awhile depending on the amount of data...\n\n');

%load file
dataraw = era_readtable('file',fullfile(pathpart,filepart));

%pull the headernames from the file
collist = dataraw.Properties.VariableNames;

%none is added so it will be presented as an option to the user (e.g., if
%there are no groups in the dataset then groups will be ignored by
%selecting none from a drop-down menu)
vcollist = collist;
vcollist{end+1} = 'none';

%era_prefs may have been defined in the function call not using the gui or
%after the chains did not converge
if ~exist('era_prefs','var')
    [era_prefs,~] = era_findprefsdata(varargin);
elseif isempty(era_prefs)
    %load default preferences for processing and viewing data
    era_prefs = era_defaults;

    %attach the current version number to era_prefs
    era_prefs.ver = era_defineversion;

    %define parameters for figure position
    era_prefs.guis.fsize = get(0,'DefaultTextFontSize');
end

%Insert the information into a data structure for holding the era data
era_data = struct;
era_data.ver = era_prefs.ver;
era_data.raw.filename = filepart;
era_data.raw.filepath = pathpart;
era_data.raw.data = dataraw;
era_data.raw.colnames = collist;
era_data.proc.collist = vcollist;

%use the loaded dateset in the gui for setting up the data to be run
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);

end

function era_startproc_gui(varargin)
%
%Inputs
% era_prefs - preferences for era toolbox
% era_data - era toolbox data structure
%
%Optional Inputs:
% procprefs - prefences for processing in stan
%  .nchains - number of chains
%  .niter - number of iterations
% inpchoices - choices for the assignment of table columns to data inputs
%  (e.g., id, measurement). Used when the preferences executes era_gui.
%
%Output
% There are no direct outputs to the Matlab workspace. The user will have
%  the option (via the gui) to go back to the era_start window or pass the
%  data to era_startproc to process the data in Stan
%

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%define parameters for figure position
figwidth = 600;
figheight = 425;

%if the assignment of which columns belong to which category has already
%been made, then use those choices. Otherwise, load the defaults and try to
%figure out the data.
if ~isfield(era_prefs.proc,'inp')
    
    %help the user set things up by checking if any generic names are present
    %in the headers that identify the columns
    era_prefs.proc.inp.id = 1;
    era_prefs.proc.inp.meas = 2;
    era_prefs.proc.inp.group = length(era_data.proc.collist);
    era_prefs.proc.inp.event = length(era_data.proc.collist);
    era_prefs.proc.inp.time = length(era_data.proc.collist);
    era_prefs.proc.inp.whichgroups = '';
    era_prefs.proc.inp.whichgroupscol = [];
    era_prefs.proc.inp.whichevents = '';
    era_prefs.proc.inp.whicheventscol = [];
    era_prefs.proc.inp.whichtimes = '';
    era_prefs.proc.inp.whichtimescol = [];

    %check for a participant header
    poss = {'Subject' 'ID' 'Participant' 'SubjID' 'Subj'};
    for i = 1:length(era_data.proc.collist)
        ind = strcmpi(era_data.proc.collist(i),poss);
        if any(ind) 
            era_prefs.proc.inp.id = i;
        end
    end

    %check for a mesurement header
    poss = {'Measurement'};
    for i = 1:length(era_data.proc.collist)
        ind = strcmpi(era_data.proc.collist(i),poss);
        if any(ind)
            era_prefs.proc.inp.meas = i;
        end
    end

    %check for a group header
    poss = {'Group'};
    for i = 1:length(era_data.proc.collist)
        ind = strcmpi(era_data.proc.collist(i),poss);
        if any(ind) 
            era_prefs.proc.inp.group = i;
        end
    end

    %check for an event type header
    poss = {'Event' 'Type'};
    for i = 1:length(era_data.proc.collist)
        ind = strcmpi(era_data.proc.collist(i),poss);
        if any(ind)
            era_prefs.proc.inp.event = i;
        end
    end

    %check for a occasion header
    poss = {'Time' 'Occasion' ' Session'};
    for i = 1:length(era_data.proc.collist)
        ind = strcmpi(era_data.proc.collist(i),poss);
        if any(ind)
            era_prefs.proc.inp.time = i;
        end
    end
end

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (figwidth/2);

%create the basic era_gui
era_gui= figure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Inputs',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' era_data.raw.filename],...
    'Position',[0 row figwidth 25]);          

%next row
row = row - rowspace*1.5;

%Print the gui headers
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Variable',...
    'Position', [figwidth/8 row figwidth/4 25]);  

uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[5*(figwidth/8) row figwidth/4 25]);

%next row
row = row - rowspace*1.5;

%print the name of the variables and place a listbox with possible options
%start with participant ID row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Participant ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(1) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.raw.colnames,'Value',era_prefs.proc.inp.id,...
    'Position', [rcol row figwidth/2 25]);  

%next row
row = row - rowspace;

%measurement row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Measurement:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(2) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.raw.colnames,'Value',era_prefs.proc.inp.meas,...
    'Position', [rcol row figwidth/2 25]); 

%next row
row = row - rowspace;
grouprow = row;

%group row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Group ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(3) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.proc.collist,'Value',era_prefs.proc.inp.group,...
    'Position', [rcol row figwidth/2.25 25]); 

uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','...',...
    'Position', [14*figwidth/15 row+2 figwidth/15 27],...
    'Callback',{@selectgroups_call,'era_prefs',era_prefs,...
    'era_data',era_data,'inplists',inplists}); 

%next row
row = row - rowspace;
eventrow = row;

%event row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Event Type:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(4) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.proc.collist,'Value',era_prefs.proc.inp.event,...
    'Position', [rcol row figwidth/2.25 25]); 

%next row
row = row - rowspace;
timerow = row;

%time row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Occasion ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(5) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.proc.collist,'Value',era_prefs.proc.inp.time,...
    'Position', [rcol row figwidth/2.25 25]); 

%Since these buttons use inplists as an input, it needed to be specified 
%after all of inputs had been placed in inplists
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','...',...
    'TooltipString','Select a subset of groups to process',...
    'Position', [14*figwidth/15 (grouprow+2) figwidth/15 27],...
    'Callback',{@selectgroups_call,'era_prefs',era_prefs,...
    'era_data',era_data,'inplists',inplists});

uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','...',...
    'TooltipString','Select a subset of events to process',...
    'Position', [14*figwidth/15 eventrow+2 figwidth/15 27],...
    'Callback',{@selectevents_call,'era_prefs',era_prefs,...
    'era_data',era_data,'inplists',inplists}); 

uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','...',...
    'TooltipString','Select a subset of occasions to process',...
    'Position', [14*figwidth/15 timerow+2 figwidth/15 27],...
    'Callback',{@selecttime_call,'era_prefs',era_prefs,...
    'era_data',era_data,'inplists',inplists}); 

%next row with extra space
row = row - rowspace*2.5;

%Create a back button that will take the user back to era_start
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 row figwidth/5 50],...
    'Callback',{@bb_call,era_gui}); 

%Create button that will check the inputs and begin processing the data
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [3*figwidth/8 row figwidth/5 50],...
    'Callback',{@era_exec,'era_prefs',era_prefs,'era_data',era_data,...
    'inplists',inplists}); 

%Create button that will display preferences
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [5*figwidth/8 row figwidth/5 50],...
    'Callback',{@era_procprefs,'era_prefs',era_prefs,'era_data',...
    era_data,'inplists',inplists}); 

%tag gui
era_gui.Tag = 'era_gui';

end


function bb_call(varargin) 
%if back button is pushed, go back to era_start

%close the era_gui
close(varargin{3});

%go back to era_start
era_start;

end

function era_procprefs(varargin)
%
%Input 
% era_prefs - toolbox preferences
% era_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The user will have
%  the option (via the gui) to go back to the era_gui or save the specified
%  preferences to be used for stan
%

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1}; 
    choices = cell2mat(get(inplists(:),'value'));
    era_prefs.proc.inp.id = choices(1);
    era_prefs.proc.inp.meas = choices(2);
    era_prefs.proc.inp.group = choices(3);
    era_prefs.proc.inp.event = choices(4);
    era_prefs.proc.inp.time = choices(5);
end

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    era_prefs.guis.pos = era_gui.Position;
    close(era_gui);
else
    era_prefs.guis.pos=[400 400 600 400];
end

%recommend that the user use at least 10k iterations if running retest
%reliability anlayses
if era_prefs.proc.inp.time ~= length(era_data.proc.collist)
    str = {'WARNING: It is recommended to run at least 10,000 iterations'};
    str{end+1} = ' if you are conducting test-retest reliability';
    str{end+1} = ' analyses. This will take a very long time, but it is ';
    str{end+1} = 'unlikely that your model will converge with fewer ';
    str{end+1} = 'iterations.';
    errordlg(str,'WARNING: Need more iterations');
end

%define space between rows and first row location
rowspace = 40;
row = era_prefs.guis.pos(4) - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (era_prefs.guis.pos(3)/2+10);

%create the basic era_prefs
era_gui= figure('unit','pix',...
  'position',era_prefs.guis.pos,...
  'menub','no',...
  'name','Specify Processing Preferences',...
  'numbertitle','off',...
  'resize','off');    

%print the gui headers
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [era_prefs.guis.pos(4)/7 row era_prefs.guis.pos(4)/3 25]);  

uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[6*(era_prefs.guis.pos(4)/7) row era_prefs.guis.pos(4)/3 25]);

%next row
row = row - rowspace*1.75;

%ask for an input for the number of stan chains
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Number of Chains:',...
    'Tooltip','Number of chains to be run in Stan (must specify at least 3)',...
    'Position', [lcol row era_prefs.guis.pos(4)/2 30]);  

newprefs.nchains = uicontrol(era_gui,'Style','edit','fontsize',era_prefs.guis.fsize,...
    'String',era_prefs.proc.nchains,...
    'Position', [rcol row+10 era_prefs.guis.pos(4)/2 30]);  

%next row
row = row - rowspace;

%input for specifying the number of iterations to run
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Number of iterations:',...
    'Tooltip','Number of iterations to run in Stan',...
    'Position', [lcol row era_prefs.guis.pos(4)/2 30]);  

newprefs.niter = uicontrol(era_gui,'Style','edit','fontsize',era_prefs.guis.fsize,...
    'String',era_prefs.proc.niter,...
    'Position', [rcol row+10 era_prefs.guis.pos(4)/2 30]);  

%next row
row = row - rowspace;

%input for specifying whether to use verbose stan output
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Verbose Stan output (print each iteration):',...
    'Tooltip','Displays Stan output in the Matlab command window',...
    'Position', [lcol row era_prefs.guis.pos(4)/2 30]);  

newprefs.verbose = uicontrol(era_gui,'Style','pop',...
    'fontsize',era_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',era_prefs.proc.verbose,...
    'Position', [rcol row+5 era_prefs.guis.pos(4)/2 30]);  

%next row
row = row - rowspace;

%input for specifying whether to view trace plots
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','View trace plots prior to saving Stan output',...
    'Tooltip','Displays trace plots to give the user the option to rerun',...
    'Position', [lcol row era_prefs.guis.pos(4)/2 30]);  

newprefs.traceplots = uicontrol(era_gui,'Style','pop',...
    'fontsize',era_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',era_prefs.proc.traceplots,...
    'Position', [rcol row+5 era_prefs.guis.pos(4)/2 30]);  

%next row with extra space
row = row - rowspace*2.5;

%Create a back button that will save inputs for preferences
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [era_prefs.guis.pos(4)/7 row era_prefs.guis.pos(4)/3 50],...
    'Callback',{@era_prefs_save,'era_prefs',era_prefs,'era_data',...
    era_data,'newprefs',newprefs}); 

%Create button that will go back to era_gui without saving
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [6*era_prefs.guis.pos(4)/7 row era_prefs.guis.pos(4)/3 50],...
    'Callback',{@era_prefs_back,'era_prefs',era_prefs,'era_data',...
    era_data}); 

%tag gui
era_gui.Tag = 'era_gui';

end

function era_prefs_save(varargin)
%if the preferences are to be saved

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find newprefs
ind = find(strcmp('newprefs',varargin),1);
newprefs = varargin{ind+1}; 

%pull the preferences and store them in era_prefs
era_prefs.proc.nchains = str2double(newprefs.nchains.String);
era_prefs.proc.niter = str2double(newprefs.niter.String);
era_prefs.proc.verbose = newprefs.verbose.Value;
era_prefs.proc.traceplots = newprefs.traceplots.Value;

%find era_gui
era_gui = findobj('Tag','era_gui');

%close era_gui
close(era_gui);

if era_prefs.proc.nchains >= 3 && era_prefs.proc.niter > 0
    %execute era_startproc_fic with the new preferences
    era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
end

%make sure the user has defined at least 3 chains and more than 0
%iterations
if era_prefs.proc.nchains < 3 
   
    errordlg('You must have at least three chains to test convergence');
    era_prefs.proc.nchains = 3;
    
    era_procprefs('era_prefs',era_prefs,'era_data',era_data);
    
end

if era_prefs.proc.niter <= 0
   
    errordlg('You must have more than zero iterations per chain');
    era_prefs.proc.niter = 1000;
    
    era_procprefs('era_prefs',era_prefs,'era_data',era_data);
    
end

end

function era_prefs_back(varargin)
%if the back button was pressed then the inputs will not be saved

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find era_gui
era_gui = findobj('Tag','era_gui');

%close era_gui
close(era_gui);

%execute era_startproc_gui with the old preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);

end

function selectgroups_call(varargin)
%select groups to be processed in the dataset
%
%Input 
% era_prefs - toolbox preferences
% era_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The inputs will be
%  stored for later processing
%

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1}; 
    choices = cell2mat(get(inplists(:),'value'));
    era_prefs.proc.inp.id = choices(1);
    era_prefs.proc.inp.meas = choices(2);
    era_prefs.proc.inp.group = choices(3);
    era_prefs.proc.inp.event = choices(4);
    era_prefs.proc.inp.time = choices(5);
end

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    era_prefs.guis.pos = era_gui.Position;
    close(era_gui);
else
    era_prefs.guis.pos=[400 400 600 400];
end

%ensure that 'none' is not selected. If it is, spit out an error and go
%back to era_startproc_gui
if strcmpi(era_data.proc.collist{era_prefs.proc.inp.group},'none')
    era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
    errordlg('No group column selected. Select a column for group',...
        'Group Column Not Defined');
    return;
end

%pull a list of groups from the file based on the input from era_gui
gnames = unique(era_data.raw.data.(era_data.proc.collist{...
    era_prefs.proc.inp.group}));

%check whether this function has been called before and whether the column
%assigned to group has changed
if isempty(era_prefs.proc.inp.whichgroupscol) || ...
        (~isempty(era_prefs.proc.inp.whichgroupscol) && ...
        era_prefs.proc.inp.whichgroupscol ~= era_prefs.proc.inp.group)
    era_prefs.proc.inp.whichgroupscol = era_prefs.proc.inp.group;
    era_prefs.proc.inp.whichgroups = gnames;
    gind = 1:length(gnames);
    
%if the column assigned to group is the same, then pull the previous inputs
elseif (~isempty(era_prefs.proc.inp.whichgroupscol) && ...
        era_prefs.proc.inp.whichgroupscol == era_prefs.proc.inp.group)   
    gind = zeros(1,length(era_prefs.proc.inp.whichgroups));
    for ii = 1:length(era_prefs.proc.inp.whichgroups)
        if ~isnumeric(gnames)
            gind(ii) = find(strcmp(era_prefs.proc.inp.whichgroups{ii},...
                gnames));
        elseif isnumeric(gnames)
            gind(ii) = find(gnames(:) == ...
                era_prefs.proc.inp.whichgroups{ii});
        end
    end
end

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

era_gui_whichgroups = figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Which Groups to Process',...
  'numbertitle','off',...
  'resize','off');  

movegui(era_gui_whichgroups,'center');

str = {'Select which groups you would like to be processed:'};

%Write text
uicontrol(era_gui_whichgroups,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);     

%bump down to next row
row = row - rowspace*8.5;

glist = uicontrol(era_gui_whichgroups,'style','list',...
     'min',0,'max',length(gnames),...
     'Value',gind,...
     'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
     'string',gnames);

%Create a button that will go back to era_proc without saving
uicontrol(era_gui_whichgroups,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_whichgroups_back_call,'era_prefs',era_prefs,...
    'era_data',era_data}); 

%Create button that will save group inputs
uicontrol(era_gui_whichgroups,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_whichgroups_save_call,'era_prefs',era_prefs,...
    'era_data',era_data,'glist',glist,'gnames',gnames});  
 
%display gui
set(era_gui_whichgroups,'Visible','on');

%tag the gui
era_gui_whichgroups.Tag = 'era_gui_whichgroups';
end

function era_whichgroups_back_call(varargin)
%go back to era_startproc and do not save which groups should be processed

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find era_gui_whichgroups
era_gui_whichgroups = findobj('Tag','era_gui_whichgroups');

%close era_gui
close(era_gui_whichgroups);

%execute era_startproc_gui with the old preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);

end

function era_whichgroups_save_call(varargin)
%go back to era_startproc and save which groups should be processed

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find glist and gnames
ind = find(strcmp('glist',varargin),1);
glist = varargin{ind+1}.Value; 
ind = find(strcmp('gnames',varargin),1);
gnames = varargin{ind+1};

groups = cell(1,length(glist));

if ~isnumeric(gnames)
    for ii = 1:length(glist)
        groups{ii} = gnames{glist(ii)};
    end   
elseif isnumeric(gnames)
    for ii = 1:length(glist)
        groups{ii} = gnames(glist(ii));
    end
end

era_prefs.proc.inp.whichgroups = groups;

%find era_gui_whichgroups
era_gui_whichgroups = findobj('Tag','era_gui_whichgroups');

%close era_gui
close(era_gui_whichgroups);

%execute era_startproc_gui with the old preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
end

function selectevents_call(varargin)
%select events to be processed in the dataset
%
%Input 
% era_prefs - toolbox preferences
% era_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The inputs will be
%  stored for later processing
%

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1}; 
    choices = cell2mat(get(inplists(:),'value'));
    era_prefs.proc.inp.id = choices(1);
    era_prefs.proc.inp.meas = choices(2);
    era_prefs.proc.inp.group = choices(3);
    era_prefs.proc.inp.event = choices(4);
    era_prefs.proc.inp.time = choices(5);
end

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    era_prefs.guis.pos = era_gui.Position;
    close(era_gui);
else
    era_prefs.guis.pos=[400 400 600 400];
end

%ensure that 'none' is not selected. If it is, spit out an error and go
%back to era_startproc_gui
if strcmpi(era_data.proc.collist{era_prefs.proc.inp.event},'none')
    era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
    errordlg('No event column selected. Select a column for event',...
        'Event Column Not Defined');
    return;
end

%pull a list of events from the file based on the input from era_gui
enames = unique(era_data.raw.data.(era_data.proc.collist{...
    era_prefs.proc.inp.event}));

%check whether this function has been called before and whether the column
%assigned to event has changed
if isempty(era_prefs.proc.inp.whicheventscol) || ...
        (~isempty(era_prefs.proc.inp.whicheventscol) && ...
        era_prefs.proc.inp.whicheventscol ~= era_prefs.proc.inp.event)
    era_prefs.proc.inp.whicheventscol = era_prefs.proc.inp.event;
    era_prefs.proc.inp.whichevents = enames;
    eind = 1:length(enames);
    
%if the column assigned to event is the same, then pull the previous inputs
elseif (~isempty(era_prefs.proc.inp.whicheventscol) && ...
        era_prefs.proc.inp.whicheventscol == era_prefs.proc.inp.event)   
    eind = zeros(1,length(era_prefs.proc.inp.whichevents));
    for ii = 1:length(era_prefs.proc.inp.whichevents)
        if ~isnumeric(enames)
            eind(ii) = find(strcmp(era_prefs.proc.inp.whichevents{ii},...
                enames));
        elseif isnumeric(enames)
            eind(ii) = find(enames(:) == ...
                era_prefs.proc.inp.whichevents{ii});
        end
    end
end

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

era_gui_whichevents = figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Which Events to Process',...
  'numbertitle','off',...
  'resize','off');  

movegui(era_gui_whichevents,'center');

str = {'Select which events you would like to be processed:'};

%Write text
uicontrol(era_gui_whichevents,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);     

%bump down to next row
row = row - rowspace*8.5;

elist = uicontrol(era_gui_whichevents,'style','list',...
     'min',0,'max',length(enames),...
     'Value',eind,...
     'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
     'string',enames);

%Create a button that will go back to era_proc without saving
uicontrol(era_gui_whichevents,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_whichevents_back_call,'era_prefs',era_prefs,...
    'era_data',era_data}); 

%Create button that will save events inputs
uicontrol(era_gui_whichevents,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_whichevents_save_call,'era_prefs',era_prefs,...
    'era_data',era_data,'elist',elist,'enames',enames});  
 
%display gui
set(era_gui_whichevents,'Visible','on');

%tag the gui
era_gui_whichevents.Tag = 'era_gui_whichevents';
end

function era_whichevents_back_call(varargin)
%go back to era_startproc and do not save which events should be processed

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find era_gui_whichevents
era_gui_whichevents = findobj('Tag','era_gui_whichevents');

%close era_gui
close(era_gui_whichevents);

%execute era_startproc_gui with the old preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);

end

function era_whichevents_save_call(varargin)
%go back to era_startproc and save which events should be processed

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find elist and enames
ind = find(strcmp('elist',varargin),1);
elist = varargin{ind+1}.Value; 
ind = find(strcmp('enames',varargin),1);
enames = varargin{ind+1};

events = cell(1,length(elist));

if ~isnumeric(enames)
    for ii = 1:length(elist)
        events{ii} = enames{elist(ii)};
    end  
elseif isnumeric(enames)
    for ii = 1:length(elist)
        events{ii} = enames(elist(ii));
    end
end

era_prefs.proc.inp.whichevents = events;

%find era_gui_whichevents
era_gui_whichevents = findobj('Tag','era_gui_whichevents');

%close era_gui
close(era_gui_whichevents);

%execute era_startproc_gui with the new preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
end

function selecttime_call(varargin)
%select occasions to be processed in the dataset
%
%Input 
% era_prefs - toolbox preferences
% era_data - toolbox data
%
%Output
% There are no direct outputs to the Matlab workspace. The inputs will be
%  stored for later processing
%

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
if ~isempty(ind)
    inplists = varargin{ind+1}; 
    choices = cell2mat(get(inplists(:),'value'));
    era_prefs.proc.inp.id = choices(1);
    era_prefs.proc.inp.meas = choices(2);
    era_prefs.proc.inp.group = choices(3);
    era_prefs.proc.inp.event = choices(4);
    era_prefs.proc.inp.time = choices(5);
end

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    era_prefs.guis.pos = era_gui.Position;
    close(era_gui);
else
    era_prefs.guis.pos=[400 400 600 400];
end

%ensure that 'none' is not selected. If it is, spit out an error and go
%back to era_startproc_gui
if strcmpi(era_data.proc.collist{era_prefs.proc.inp.time},'none')
    era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
    errordlg('No occasion column selected. Select a column for occasion',...
        'Occasion Column Not Defined');
    return;
end

%pull a list of occasions from the file based on the input from era_gui
tnames = unique(era_data.raw.data.(era_data.proc.collist{...
    era_prefs.proc.inp.time}));

%check whether this function has been called before and whether the column
%assigned to event has changed
if isempty(era_prefs.proc.inp.whichtimescol) || ...
        (~isempty(era_prefs.proc.inp.whichtimescol) && ...
        era_prefs.proc.inp.whichtimescol ~= era_prefs.proc.inp.time)
    era_prefs.proc.inp.whichtimescol = era_prefs.proc.inp.time;
    era_prefs.proc.inp.whichtimes = tnames;
    tind = 1:length(tnames);
    
%if the column assigned to occasion is the same, then pull the previous inputs
elseif (~isempty(era_prefs.proc.inp.whichtimescol) && ...
        era_prefs.proc.inp.whichtimescol == era_prefs.proc.inp.time)   
    tind = zeros(1,length(era_prefs.proc.inp.whichtimes));
    for ii = 1:length(era_prefs.proc.inp.whichtimes)
        if ~isnumeric(tnames)
            tind(ii) = find(strcmp(era_prefs.proc.inp.whichtimes{ii},...
                tnames));
        elseif isnumeric(tnames)
            tind(ii) = find(tnames(:) == ...
                era_prefs.proc.inp.whichtimes{ii});
        end
    end
end

%define parameters for figure position
figwidth = 500;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

era_gui_whichtimes = figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Which Occasions to Process',...
  'numbertitle','off',...
  'resize','off');  

movegui(era_gui_whichtimes,'center');

str = {'Select which events you would like to be processed:'};

%Write text
uicontrol(era_gui_whichtimes,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 25]);     

%bump down to next row
row = row - rowspace*8.5;

tlist = uicontrol(era_gui_whichtimes,'style','list',...
     'min',0,'max',length(tnames),...
     'Value',tind,...
     'Position',[.5*figwidth/4 row 3*figwidth/4 figheight/2],...
     'string',tnames);

%Create a button that will go back to era_proc without saving
uicontrol(era_gui_whichtimes,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_whichtimes_back_call,'era_prefs',era_prefs,...
    'era_data',era_data}); 

%Create button that will save events inputs
uicontrol(era_gui_whichtimes,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [4.5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_whichtimes_save_call,'era_prefs',era_prefs,...
    'era_data',era_data,'tlist',tlist,'tnames',tnames});  
 
%display gui
set(era_gui_whichtimes,'Visible','on');

%tag the gui
era_gui_whichtimes.Tag = 'era_gui_whichtimes';
end

function era_whichtimes_back_call(varargin)
%go back to era_startproc and do not save which occasions should be processed

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find era_gui_whichevents
era_gui_whichtimes = findobj('Tag','era_gui_whichtimes');

%close era_gui
close(era_gui_whichtimes);

%execute era_startproc_gui with the old preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);

end

function era_whichtimes_save_call(varargin)
%go back to era_startproc and save which events should be processed

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find tlist and tnames
ind = find(strcmp('tlist',varargin),1);
tlist = varargin{ind+1}.Value; 
ind = find(strcmp('tnames',varargin),1);
tnames = varargin{ind+1};

times = cell(1,length(tlist));

if ~isnumeric(tnames)
    for ii = 1:length(tlist)
        times{ii} = tnames{tlist(ii)};
    end  
elseif isnumeric(tnames)
    for ii = 1:length(tlist)
        times{ii} = tnames(tlist(ii));
    end
end

era_prefs.proc.inp.whichtimes = times;

%find era_gui_whichevents
era_gui_whichtimes = findobj('Tag','era_gui_whichtimes');

%close era_gui
close(era_gui_whichtimes);

%execute era_startproc_gui with the new preferences
era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
end

function era_exec(varargin)
%if execute button is pushed, parse input to run in era_relwrap
%
%Input
% era_prefs - toolbox preferences
% era_data - toolbox data
% inplists - list of input choices from era_startproc_gui
%
%Output
% No variables outputted to the Matlab workspace
% A directory will be created where the temporary files will be saved for
%  running the Stan model (stan creates various temp files). After the
%  model is done running (executed with era_computerel), the temporary
%  directory and files will be removed.
% The user will also be asked to specify the location of where the .mat
%  file that contains the Stan results should be saved. This .mat file will
%  be used by era_startview, which will be automatically pulled up after
%  stan has finished.

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inplists
ind = find(strcmp('inplists',varargin),1);
inplists = varargin{ind+1}; 
choices = cell2mat(get(inplists(:),'value'));
era_prefs.proc.inp.id = choices(1);
era_prefs.proc.inp.meas = choices(2);
era_prefs.proc.inp.group = choices(3);
era_prefs.proc.inp.event = choices(4);
era_prefs.proc.inp.time = choices(5);

%parse inputs
era_data.proc.idheader = char(era_data.proc.collist(choices(1)));
era_data.proc.measheader = char(era_data.proc.collist(choices(2)));
    
if choices(3) ~= length(era_data.proc.collist)
    era_data.proc.groupheader = char(era_data.proc.collist(choices(3)));
elseif choices(3) == length(era_data.proc.collist)
    era_data.proc.groupheader = '';
end

if choices(4) ~= length(era_data.proc.collist)
    era_data.proc.eventheader = char(era_data.proc.collist(choices(4)));
elseif choices(4) == length(era_data.proc.collist)
    era_data.proc.eventheader = '';
end

if choices(5) ~= length(era_data.proc.collist)
    era_data.proc.timeheader = char(era_data.proc.collist(choices(5)));
elseif choices(5) == length(era_data.proc.collist)
    era_data.proc.timeheader = '';
end

%check whether particular events were specified to process. If not, process
%all event types. Also, if event was changed, overwrite the old and replace
%with all the new event types.
if (isempty(era_prefs.proc.inp.whicheventscol) && ...
        ~strcmpi(era_data.proc.collist{era_prefs.proc.inp.event},'none')) || ...
        (~isempty(era_prefs.proc.inp.whicheventscol) && ...
        era_prefs.proc.inp.whicheventscol ~= era_prefs.proc.inp.event)
    
    %pull a list of events from the file based on the input from era_gui
    enames = unique(era_data.raw.data.(era_data.proc.collist{...
        era_prefs.proc.inp.event}));

    era_prefs.proc.inp.whicheventscol = era_prefs.proc.inp.event;
    era_prefs.proc.inp.whichevents = enames;   
    
end

%check whether particular group were specified to process. If not, process
%all group types. Also, if group was changed, overwrite the old and replace
%with all the new group types.
if (isempty(era_prefs.proc.inp.whichgroupscol) &&...
        ~strcmpi(era_data.proc.collist{era_prefs.proc.inp.group},'none')) || ...
        (~isempty(era_prefs.proc.inp.whichgroupscol) && ...
        era_prefs.proc.inp.whichgroupscol ~= era_prefs.proc.inp.group)
    
    %pull a list of groups from the file based on the input from era_gui
    gnames = unique(era_data.raw.data.(era_data.proc.collist{...
        era_prefs.proc.inp.group}));
    
    era_prefs.proc.inp.whichgroupscol = era_prefs.proc.inp.group;
    era_prefs.proc.inp.whichgroups = gnames;  
    
end

%check whether particular occasions were specified to process. If not, process
%all occasions. Also, if occasion was changed, overwrite the old and replace
%with all the new occasion types.
if (isempty(era_prefs.proc.inp.whichtimescol) &&...
        ~strcmpi(era_data.proc.collist{era_prefs.proc.inp.time},'none')) || ...
        (~isempty(era_prefs.proc.inp.whichtimescol) && ...
        era_prefs.proc.inp.whichtimescol ~= era_prefs.proc.inp.time)
    
    %pull a list of occasions from the file based on the input from era_gui
    tnames = unique(era_data.raw.data.(era_data.proc.collist{...
        era_prefs.proc.inp.time}));
    
    era_prefs.proc.inp.whichtimescol = era_prefs.proc.inp.time;
    era_prefs.proc.inp.whichtimes = tnames;  
    
end

%put the information about which events, groups, and occasions to process 
%in era_data
if ~strcmpi(era_data.proc.collist{era_prefs.proc.inp.event},'none')
    era_data.proc.whichevents = era_prefs.proc.inp.whichevents;
else
    era_data.proc.whichevents = '';
end
if ~strcmpi(era_data.proc.collist{era_prefs.proc.inp.group},'none')
    era_data.proc.whichgroups = era_prefs.proc.inp.whichgroups;
else
    era_data.proc.whichgroups = '';
end
if ~strcmpi(era_data.proc.collist{era_prefs.proc.inp.time},'none')
    era_data.proc.whichtimes = era_prefs.proc.inp.whichtimes;
else
    era_data.proc.whichtimes = '';
end

%find era_gui
era_gui = findobj('Tag','era_gui');

%close era_gui
close(era_gui);

%check if there are duplicates (other than 'none') (e.g., group and event
%do not refer to the same column in the data)
[n, bin] = histc(choices, unique(choices));
multiple = find(n > 1);
ind = find(ismember(bin, multiple));
if ~isempty(ind)
    probcol = {};
    for i = 1:length(ind)
        if choices(ind(i)) ~= length(era_data.proc.collist)
            probcol(end+1) = era_data.proc.collist(ind(i)); %#ok<AGROW>
        end
    end
    if ~isempty(probcol)
        dlg = {'Duplicate variable names were not provided for '};
        for i = 1:length(probcol)
            dlg{end+1} = probcol{i}; %#ok<AGROW>
        end
        dlg{end+1} = 'When selecting column headers, please select unique names';
        errordlg(dlg, 'Unique variable names not provided');
        
        %take the user back to era_startproc_gui
        era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
        return;
    end
end

%make sure the data in the column labeled measheader are actually numeric.
%If the data are not numeric, give an error and take the user back to
%era_startproc_gui
if ~isnumeric(era_data.raw.data.(era_data.proc.measheader))
    dlg = {'The data column specified by';...
        era_data.proc.measheader; 'does not contain numeric data'; ...
        'Only numeric data are allowed in the Measurement variable';...
        'Please select a different column in the dataset for Measurement'};
    errordlg(dlg, 'Measurement data not numeric');

    %take the user back to era_startproc_gui
    era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
    return;
end

%cmdstan cannot handle paths with white space. User will be required to
%provide a path to working directory that does not include a space.
%space will be set to 1 and will be changed if the path does not include
%whitespace
space = 1;
while space == 1
    
    %prompt the user to indicate where the output from stan should be saved
    [era_data.proc.savename, era_data.proc.savepath] = uiputfile(...
        fullfile(era_data.raw.filepath,'*.erat'),...
        'Save output file as');
    
    if any(isspace(era_data.proc.savename)) ||...
            any(isspace(era_data.proc.savepath))
       str = {};
       str{end+1} = 'Cmdstan cannot handle file paths or filesnames with whitespace';
       str{end+1} = 'Please specify a file path and filename without whitespace';
        errordlg(str,'Whitespace detected');
        
    elseif ~any(isspace(era_data.proc.savename)) &&...
            ~any(isspace(era_data.proc.savepath))
        
        space = 0;
        
    end
    
    %if the user does not select a file, then take the user back to 
    %era_startproc_gui    
    if era_data.proc.savename == 0 
        errordlg('No file selected','File Error');
        era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
    end
    
end %while space == 1

%in the event that the user cancels and does not specify and path and
%filename
if isempty(era_data.proc.savename) || isempty(era_data.proc.savepath)
    %take the user back to era_startproc_gui
        era_startproc_gui('era_prefs',era_prefs,'era_data',era_data);
end

%now that the headers have been specified, set up a data table to be passed
%to era_computerel for analysis
era_data.proc.data = era_loadfile('era_prefs',era_prefs,...
    'era_data',era_data);

%pass the data to be setup for processing
era_data = era_computerelwrap('era_prefs',era_prefs,'era_data',era_data);

%only continue to view the data if the reliability data were put into the
%era_data structure. If the chains did not converge there will be no rel
%field
if isfield(era_data,'rel')
    %take the user to era_startview for viewing the processed data
    era_startview('era_prefs',era_prefs,'era_data',era_data);
end

end