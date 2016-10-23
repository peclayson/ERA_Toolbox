function era_startproc(varargin)
%
%Initiate Matlab gui to begin processing data in Stan
%
%
%Last Updated 10/22/16
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

%check if era_gui is open. If the user executes era_startproc and skips
%era_start then there will be no gui to close.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    %Grab preferences if they exist
    if ~isempty(varargin)

        %check if for era_prefs
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
    errordlg('No file selected','File Error');
    era_start;
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
figheight = 400;

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

    %check for a event type header
    poss = {'Event' 'Type'};
    for i = 1:length(era_data.proc.collist)
        ind = strcmpi(era_data.proc.collist(i),poss);
        if any(ind)
            era_prefs.proc.inp.event = i;
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

%group row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Group ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(3) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.proc.collist,'Value',era_prefs.proc.inp.group,...
    'Position', [rcol row figwidth/2 25]); 

%next row
row = row - rowspace;

%event row
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Event Type:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(4) = uicontrol(era_gui,'Style','pop','fontsize',era_prefs.guis.fsize,...
    'String',era_data.proc.collist,'Value',era_prefs.proc.inp.event,...
    'Position', [rcol row figwidth/2 25]); 

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
end

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    era_prefs.guis.pos = era_gui.Position;
    close(era_gui);
else
    era_prefs.guis.pos=[400 400 600 400];
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

%input for number of stan iterations
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

%input for number of stan iterations
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Verbose Stan output (print each iteration):',...
    'Tooltip','Displays Stan output in the Matlab command window',...
    'Position', [lcol row era_prefs.guis.pos(4)/2 30]);  

newprefs.verbose = uicontrol(era_gui,'Style','pop',...
    'fontsize',era_prefs.guis.fsize,'String',{'No';'Yes'},...
    'Value',era_prefs.proc.verbose,...
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
            probcol(end+1) = era_data.proc.collist(ind(i));
        end
    end
    if ~isempty(probcol)
        dlg = {'Duplicate variable names were not provided for '};
        for i = 1:length(probcol)
            dlg{end+1} = probcol{i};
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
        'Where would you like to save the output files?');
    
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

%take the user to era_startview for viewing the processed data
era_startview('era_prefs',era_prefs,'era_data',era_data);

end