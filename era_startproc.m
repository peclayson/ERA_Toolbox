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

%check if era_gui is open. If the user executes era_startproc and skips
%era_start then there will be no gui to close.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

%ask the user to identify the data file to be loaded
[filepart, pathpart] = uigetfile({'*.xlsx','Excel File (.xlsx)';'*.xls',...
        'Excel File 97-2003 (.xls)';'*.csv',...
        'Comma-Separated Vale File (.csv)';'*.dat',...
        'Tab-Delimited Text (.dat)';'*.txt',...
        'Raw Text File (.txt)';'*.ods',...
        'OpenDocument Spreadsheet (.ods)'},'Data');

%if the user does not select a file, then take the user back to era_start    
if filepart == 0 
    errordlg('No file selected','File Error');
    era_start;
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
collist{end+1} = 'none';

%use the loaded dateset in the gui for setting up the data to be run
era_startproc_fig(collist,filepart,pathpart,dataraw);

end

function era_startproc_fig(collist,filepart,pathpart,dataraw,varargin)
%
%Input
% collist - a list of the headers of the loaded dataset
% filepart - the name of the dataset
% pathpart - the path where the dataset is located
% dataraw - the loaded dataset
% procprefs - processing preferences to be passed to Stan
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

%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value\n'));
    end
 
    %check if procprefs has already been defined
    ind = find(strcmp('procprefs',varargin),1);
    if ~isempty(ind)
        procprefs = varargin{ind+1}; 
    else
        %define default preferences for processing data in stan
        procprefs.nchains = 5;
        procprefs.niter = 250;
    end

    %check if inpchoices has already been defined
    ind = find(strcmp('inpchoices',varargin),1);
    if ~isempty(ind)
        inpchoices = varargin{ind+1}; 
    else
        inpchoices = [];
    end

elseif isempty(varargin) %just in case none of the optional inputs are used
    
    procprefs.nchains = 3;
    procprefs.niter = 500;
    inpchoices = [];
    
end %if ~isempty(varargin)

%define parameters for figure position
figwidth = 600;
figheight = 400;
collist_nonone = collist;
collist_nonone(end) = [];

%if the assignment of which columns belong to which category has already
%been made, then use those choices. Otherwise, load the defaults and try to
%figure out the data.
if ~isempty(inpchoices)
    
    defls.id = inpchoices(1);
    defls.meas = inpchoices(2);
    defls.group = inpchoices(3);
    defls.event = inpchoices(4);
    
elseif isempty(inpchoices)

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
    'Position', [figwidth/8 row figwidth/5 50],...
    'Callback',{@bb_call,era_gui}); 

%Create button that will check the inputs and begin processing the data
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [3*figwidth/8 row figwidth/5 50],...
    'Callback',{@era_exec,inplists,collist,filepart,pathpart,...
    dataraw,procprefs,era_gui}); 

%Create button that will display preferences
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [5*figwidth/8 row figwidth/5 50],...
    'Callback',{@era_procprefs,inplists,collist,filepart,pathpart,...
    dataraw,era_gui,procprefs}); 

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
%Input (varargin #)
% inpchoices(3) - a list of the choices made for assigning headers to
%  categories (id, measurement, etc)
% collist(4) - a list of the headers of the loaded dataset
% filepart(5) - the name of the dataset
% pathpart(6) - the path where the dataset is located
% dataraw(7) - the loaded dataset
% era_gui(8) - handle for era_gui
% procprefs(9) - processing preferences to be passed to Stan
%
%Output
% There are no direct outputs to the Matlab workspace. The user will have
%  the option (via the gui) to go back to the era_gui or save the specified
%  preferences to be used for stan
%

%parse inputs
if any(ishandle(varargin{3}))
    h_era_gui.inpchoices = cell2mat(get(varargin{3},'value'));
else
    h_era_gui.inpchoices = varargin{3};
end

h_era_gui.collist = varargin{4};
h_era_gui.filepart = varargin{5};
h_era_gui.pathpart = varargin{6};
h_era_gui.dataraw = varargin{7};
initialprefs = varargin{9};
newprefs = varargin{9};

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    pos = varargin{8}.Position;
    close(era_gui);
else
    pos=[400 400 600 400];
end

%define space between rows and first row location
rowspace = 35;
row = pos(4) - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (pos(3)/2+10);

%create the basic era_prefs
era_prefs= figure('unit','pix',...
  'position',pos,...
  'menub','no',...
  'name','Specify Processing Preferences',...
  'numbertitle','off',...
  'resize','off');    

%print the gui headers
uicontrol(era_prefs,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [pos(4)/7 row pos(4)/3 25]);  

uicontrol(era_prefs,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[6*(pos(4)/7) row pos(4)/3 25]);

%next row
row = row - rowspace*1.75;

%ask for an input for the number of stan chains
uicontrol(era_prefs,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Number of Chains:',...
    'Position', [lcol row pos(4)/2 35]);  

newprefs.nchains = uicontrol(era_prefs,'Style','edit','fontsize',14,...
    'String',newprefs.nchains,...
    'Position', [rcol row pos(4)/2 25]);  

%next row
row = row - rowspace;

%input for number of stan iterations
uicontrol(era_prefs,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Number of iterations:',...
    'Position', [lcol row pos(4)/2 25]);  

newprefs.niter = uicontrol(era_prefs,'Style','edit','fontsize',14,...
    'String',newprefs.niter,...
    'Position', [rcol row pos(4)/2 25]);  


%next row with extra space
row = row - rowspace*2.5;

%Create a back button that will save inputs for preferences
uicontrol(era_prefs,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [pos(4)/7 row pos(4)/3 50],...
    'Callback',{@era_prefs_save,era_prefs,newprefs,h_era_gui}); 

%Create button that will go back to era_gui without saving
uicontrol(era_prefs,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [6*pos(4)/7 row pos(4)/3 50],...
    'Callback',{@era_prefs_back,era_prefs,initialprefs,h_era_gui}); 

end

function era_prefs_save(varargin)
%if the preferences are to be saved

%pull the preferences from the gui
procprefs.nchains = str2double(varargin{4}.nchains.String);
procprefs.niter = str2double(varargin{4}.niter.String);

%close era_prefs gui
close(varargin{3});

%pull the era_gui data
h_era_gui = varargin{5};

if procprefs.nchains >= 3 && procprefs.niter > 0
    %execute era_startproc_fic with the new preferences
    era_startproc_fig(h_era_gui.collist,h_era_gui.filepart,...
        h_era_gui.pathpart,h_era_gui.dataraw,'procprefs',procprefs,...
        'inpchoices',h_era_gui.inpchoices);
end

%make sure the user has defined at least 3 chains and more than 0
%iterations
if procprefs.nchains < 3 
   
    errordlg('You must have at least three chains to test convergence');
    procprefs.nchains = 3;
    
    era_procprefs([],[],h_era_gui.inpchoices,h_era_gui.collist,...
        h_era_gui.filepart,h_era_gui.pathpart,...
        h_era_gui.dataraw,[],procprefs);
    
end

if procprefs.niter <= 0
   
    errordlg('You must have more than zero iterations per chain');
    procprefs.niter = 500;
    
    era_procprefs([],[],h_era_gui.inpchoices,h_era_gui.collist,...
        h_era_gui.filepart,h_era_gui.pathpart,...
        h_era_gui.dataraw,[],procprefs);
    
end

end

function era_prefs_back(varargin)
%if the back button was pressed then the inputs will not be saved

%pull the old preferences
procprefs.nchains = varargin{4}.nchains;
procprefs.niter = varargin{4}.niter;

%close era_prefs gui
close(varargin{3});

%pull the era_gui data
h_era_gui = varargin{5};

%execute era_startproc_fig with the old preferences
era_startproc_fig(h_era_gui.collist,h_era_gui.filepart,...
    h_era_gui.pathpart,h_era_gui.dataraw,'procprefs',procprefs,...
    'inpchoices',h_era_gui.inpchoices);

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
procprefs = varargin{8};

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
close(varargin{9});

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
if ~exist(fullfile(savepath,'Temp_StanFiles'), 'dir')
  mkdir(savepath,'Temp_StanFiles');
end
origdir = cd(fullfile(savepath,'Temp_StanFiles'));

%set initial state of rerun to 1
%this will run era_computerel
%if chains properly converged then rerun will be changed to 0 and the while
%loop will be exited
rerun = 1;

%whether chains converged will be checked each time era_computerel is run
while rerun ~= 0

    %pass the data to era_computerel for analysis
    REL = era_computerel('data',dataout,'chains',procprefs.nchains,...
        'iter',procprefs.niter);
    
    %check convergence of chains
    RELout = era_checkconv(REL);
    
    if RELout.out.conv.converged == 0
       
        era_reruncheck;
        era_gui = findobj('Tag','era_gui');
        rerun = guidata(era_gui);
        close(era_gui);
        
    else 
        %if chains converged, do no rerun.
        rerun = 0;
        
    end
    
    %if convergence was not met and the user would like to rerun the model,
    %double the number of iterations (if the doubled number is less than
    %1000, then run 1000 iterations)
    if rerun == 1
        procprefs.niter = procprefs.niter * 2;
        if procprefs.niter < 1000
            procprefs.niter = 1000;
        end
        fprintf('\nIncreasing number of iterations to %d\n',...
            procprefs.niter);
    end
end

%sometimes the era_gui doesn't close
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end
    
%change the working directory back to the original directory 
cd(origdir);

%attempt to remove the temporary directory and its files
try
    
    rmdir(fullfile(savepath,'Temp_StanFiles'),'s');

catch
    
    fprintf('\n\nTemporary directory could not be removed.');
    fprintf('\nPath: %s\n',savepath);
end

%if the chains did not converge send the user back to era_startproc_fig
if RELout.out.conv.converged == 0

    era_startproc_fig(collist,filepart,pathpart,dataraw,'procprefs',...
        procprefs,'inpchoices',choices);
    return
    
else
    
    %chains converged. Let the user know.
    %str = 'Models successfully converged with %d chains and %d iterations';
    %fprintf(strcat('\n',str,'\n'),procpref.nchains,procprefs.niter);
    
end

fprintf('\nSaving Processed Data...\n\n');

save(fullfile(savepath,savename),'RELout');

%take the user to era_startview for viewing the processed data
era_startview('file',fullfile(savepath,savename));

end


