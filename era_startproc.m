function era_startproc(varargin)
%
%Initiate Matlab gui to begin processing data in Stan
%
%
%Last Updated 3/5/16
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
%  This .mat will automatically be loaded into era_startview for viewing the
%  data
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
% by Peter Clayson (3/5/16)
% peter.clayson@gmail.com
%

%close the era_start gui
close(varargin{3});

%ask the user to identify the data file to be loaded
[filepart, pathpart] = uigetfile({'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},'Data');

%if the user does not select a file, then take the user back to era_start    
if filepart == 0 
    errordlg('No file selected','File Error');
    era_start;
end

fprintf('\n\nLoading Data...');
fprintf('\nThis may take awhile depending on the amount of data...\n\n');
%load file
dataraw = readtable(fullfile(pathpart,filepart));

%pull the headernames from the file
collist = dataraw.Properties.VariableNames;

%none is added so it will be presented as an option to the user (e.g., if
%there are no groups in the dataset then groups will be ignored by
%selecting none from a drop-down menu)
collist{end+1} = 'none';

%use the loaded dateset in the gui for setting up the data to be run
era_startproc_fig(collist,filepart,pathpart,dataraw);

end

function era_startproc_fig(collist,filepart,pathpart,dataraw)
%
%Input
% collist - a list of the headers of the loaded dataset
% filepart - the name of the dataset
% pathpart - the path where the dataset is located
% dataraw - the loaded dataset
%
%Output
% There are no direct outputs to the Matlab workspace. The user will have
%  the option (via the gui) to go back to the era_start window or pass the
%  data to era_startproc to process the data in Stan
%

%define parameters for figure position
figwidth = 550;
figheight = 400;
collist_nonone = collist;
collist_nonone(end) = [];

%help the user set things up by checking if any generic names are present
%in the headers that identify the columns
defls = struct();
defls.id = 1;
defls.meas = 2;
defls.group = length(collist);
defls.event = length(collist);

%check for a participant header
poss = {'Subject' 'ID' 'Participant' 'SubjID' 'Subj'};
for i = 1:length(collist_nonone)
    ind = strcmpi(collist_nonone(i),poss);
    if sum(ind) == 1
        defls.id = i;
    end
end

%check for a mesurement header
poss = {'Measurement'};
for i = 1:length(collist_nonone)
    ind = strcmpi(collist_nonone(i),poss);
    if sum(ind) == 1
        defls.meas = i;
    end
end

%check for a group header
poss = {'Group'};
for i = 1:length(collist_nonone)
    ind = strcmpi(collist_nonone(i),poss);
    if sum(ind) == 1
        defls.group = i;
    end
end

%check for a group header
poss = {'Event' 'Type'};
for i = 1:length(collist_nonone)
    ind = strcmpi(collist_nonone(i),poss);
    if sum(ind) == 1
        defls.event = i;
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
uicontrol(era_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' filepart],...
    'Position',[0 row figwidth 25]);          

%next row
row = row - rowspace*1.5;

%Print the gui headers
uicontrol(era_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Variable',...
    'Position', [figwidth/8 row figwidth/4 25]);  

uicontrol(era_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[5*(figwidth/8) row figwidth/4 25]);

%next row
row = row - rowspace*1.5;

%print the name of the variables and place a listbox with possible options
%start with participant ID row
uicontrol(era_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Participant ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(1) = uicontrol(era_gui,'Style','pop','fontsize',14,...
    'String',collist_nonone,'Value',defls.id,...
    'Position', [rcol row figwidth/2 25]);  

%next row
row = row - rowspace;

%measurement row
uicontrol(era_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Measurement:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(2) = uicontrol(era_gui,'Style','pop','fontsize',14,...
    'String',collist_nonone,'Value',defls.meas,...
    'Position', [rcol row figwidth/2 25]); 

%next row
row = row - rowspace;

%group row
uicontrol(era_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Group ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(3) = uicontrol(era_gui,'Style','pop','fontsize',14,...
    'String',collist,'Value',defls.group,...
    'Position', [rcol row figwidth/2 25]); 

%next row
row = row - rowspace;

%event row
uicontrol(era_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Event Type:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(4) = uicontrol(era_gui,'Style','pop','fontsize',14,...
    'String',collist,'Value',defls.event,...
    'Position', [rcol row figwidth/2 25]); 

%next row with extra space
row = row - rowspace*2.5;

%Create a back button that will take the user back to era_start
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 row figwidth/4 50],...
    'Callback',{@bb_call,era_gui}); 

%Create button that will check the inputs and begin processing the data
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [5*figwidth/8 row figwidth/4 50],...
    'Callback',{@era_exec,inplists,collist,filepart,pathpart,dataraw,era_gui}); 


end


function bb_call(varargin) 
%if back button is pushed, go back to era_start

%close the era_gui
close(varargin{3});

%go back to era_start
era_start;
   
end


function era_exec(varargin)
%if execute button is pushed, analyze the loaded data
%
%Input
% varargin containing
%  3-inputs specifying the columns of data to be analyzed
%  4-the list of possible column headers that were chosen from
%  5-the name of the file to be analyzed
%  6-the path to dataset
%  7-the loaded dataset
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

%parse inputs
inp = varargin{3};
collist = varargin{4};
filepart = varargin{5};
pathpart = varargin{6};
dataraw = varargin{7};

choices = cell2mat(get(inp(:),'value'));

%check if there are duplicates (other than 'none') (e.g., group and event
%do not refer to the same column in the data)
[n, bin] = histc(choices, unique(choices));
multiple = find(n > 1);
ind = find(ismember(bin, multiple));
if ~isempty(ind)
    probcol = {};
    for i = 1:length(ind)
        if choices(ind(i)) ~= length(collist)
            probcol(end+1) = collist(ind(i));
        end
    end
    if ~isempty(probcol)
        dlg = {'Duplicate variable names were not provided for ';...
            probcol; ...
            'When selecting column headers, please select unique names'};
        errordlg(dlg, 'Unique variable names not provided');
        
        %take the user back to era_startproc_fig
        era_startproc_fig(collist,filepart,pathpart,dataraw);
    end
end

%parse inputs
idheader = char(collist(choices(1)));
measheader = char(collist(choices(2)));

if choices(3) ~= length(collist)
    groupheader = char(collist(choices(3)));
elseif choices(3) == length(collist)
    groupheader = '';
end

if choices(4) ~= length(collist)
    eventheader = char(collist(choices(4)));
elseif choices(4) == length(collist)
    eventheader = '';
end

%close the gui
close(varargin{8});

%prompt the user to indicate where the output from stan should be saved
[savename, savepath] = uiputfile(fullfile(pathpart,'*.mat'),...
    'Where would you like to save the output files?');

%if the user does not select a file, then take the user back to 
%era_startproc_fig    
if filepart == 0 
    errordlg('No file selected','File Error');
    era_startproc_fig(collist,filepart,pathpart,dataraw);
end

%now that the headers have been specified, set up a data table to be passed
%to era_computerel for analysis
dataout = era_loadfile('file',fullfile(pathpart,filepart),...
    'idcol',idheader,'meascol',measheader,'groupcol',groupheader,...
    'eventcol',eventheader,'dataraw',dataraw);

%Change working dir for temporary Stan files
mkdir(savepath,'Temp_StanFiles');
origdir = cd(fullfile(savepath,'Temp_StanFiles'));

%pass the data to era_computerel for analysis
RELout = era_computerel('data',dataout);

%change the working directory back to the original directory 
cd(origdir);

%attempt to remove the temporary directory and its files
try
    
    rmdir(fullfile(savepath,'Temp_StanFiles'),'s');

catch
    
    fprintf('\n\nTemporary directory could not be removed.');
    sprintf('\nPath: %s\n',savepath);
end

fprintf('\nSaving Processed Data...\n\n');

save(fullfile(savepath,savename),'RELout');

%take the user to era_startview for viewing the processed data
era_startview('file',fullfile(savepath,savename));

end


