function era_startview(varargin)
%Prepares the data for viewing and lets user specify which tables and
%figures to present
%
%era_startview('file','/Users/REL/SomeERAData.mat')
%
%Last Updated 7/24/16
%
%Required Inputs:
% No inputs are required.
%
%Optional Inputs:
% file - file of data processed using era_computerel. This optional input
%  is used by era_startproc to provide an easy transition from processing
%  to viewing without the user having to re-select a file.
%
%Output:
% No data are outputted to the Matlab command window. However, the user
%  will have the option of saving various figures and plots that will be 
%  created by era_relfigures, which is executed by this gui 

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
% changes consistent with ERA Toolbox file format (extension: .erat)
%
%4/27/16 PC
% add check to ensure that dependability estimate provided by user is
%  numeric and between 0 and 1
%
%7/20/16 PC
% consolidate option for requesting tables for ICCs and stddevs
%
%7/21/16 PC
% changes associated with adding era_prefs and era_data
%
%7/24/16 PC
% use era_data as input into era_relfigures

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%see if the file for the figures and tables has been specified in
%varargin
if ~isempty(varargin) && (isempty(era_data) && isempty(era_prefs))
    
    %check if data file has been provided
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = varargin{ind+1};
        [pathpart,filepart] = fileparts(file);
    end

end %if ~isempty(varargin)

%check if era_gui is open. If the user executes era_startproc and skips
%era_start then there will be no gui to close.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

%if the file was not specified, prompt the user to indicate where the file
%is located.
if ~exist('file','var') && isempty(era_data)
    [filepart, pathpart] = uigetfile({'*.erat',...
        'ERA Toolbox files (*.erat)'},'Data');

    if filepart == 0 
        errordlg('No file selected','File Error');
        era_start;
        return;
    end

    fprintf('\n\nLoading Data...\n\n');
    
    %load data
    load(fullfile(pathpart,filepart),'-mat');
   
end

%if era_prefs does not exist, load the default preferences. If this window
%was not gotten to using era_start, era_prefs will need to be defined
if isempty(era_prefs)
    era_prefs = era_defaults;
    era_prefs.ver = era_defineversion;
end

%create a gui to allow the user to specify what aspects of the data will be
%viewed
era_startview_fig('era_prefs',era_prefs,'era_data',era_data);

end

function era_startview_fig(varargin)
%Input
% era_data - ERA Toolbox data structure array
% era_prefs - ERA Toolbox preferences structure array
%
%Output
% No variables will be outputted to the Matlab workspace. Based on the
%  inputs from this gui, era_relfigures will be executed to display various
%  figures and tables (for more information about the tables and figures
%  see the user manual for the ERA toolbox)

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

%define parameters for figure position
figwidth = 550;
figheight = 550;

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

%create the gui
era_gui= figure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Inputs',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' era_data.rel.filename],...
    'Position',[0 row figwidth 25]);          

%next row
row = row - (rowspace*1.5);

%Print the text for dependability cutoff with a box for the user to specify
%the input
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Dependability Cutoff:',...
    'Position', [lcol row figwidth/4 25]);  

inputs.h(1) = uicontrol(era_gui,'Style','edit','fontsize',era_prefs.guis.fsize,...
    'String',era_prefs.view.depvalue,...
    'Position', [rcol+5 row+6 figwidth/4 25]);  

%next row
row = row - rowspace-10;

%indicate that a checked box means yes
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Checked = YES',...
    'Position', [rcol+5 row figwidth/4 25]);  

%next row
row = row - rowspace;

%increase distance between rows as some descriptions take up more than one
%line
rowspace = 50;
rcol = (figwidth/4)*3;

%dependability with increasing trials
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Would you like to plot Number of Trials v Dependability?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(2) = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.plotdep,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%plot ICCs
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like to plot intraclass correlation coefficients?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(3) = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.ploticc,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%dependability cutoff table
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like a table of event specific dependability information?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(4) = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.inctrltable,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%overall dependability table
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like a table of overall dependability information?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(5) = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.overalltable,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%between- and within-person standard deviation tables
str = ['Would you like a table of the overall relative '...
    'sizes of sources of variance and ICCs?'];
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(6) = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.showstddevt,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row
row = row - rowspace;

%plot between-person standard deviations
str = ['Would you like a figure showing the between-person standard '...
    'deviations?'];
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(7) = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.showstddevf,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row for buttons
row = row - rowspace*1.5;

%Create a back button that will take the user back to era_start
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 row figwidth/5 50],...
    'Callback',{@era_svb,era_gui}); 

%Create button that will check the inputs and begin processing the data
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [3*figwidth/8 row figwidth/5 50],...
    'Callback',{@era_svh,'era_prefs',era_prefs,'era_data',era_data,...
    'inputs',inputs}); 

%Create button that will display preferences
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [5*figwidth/8 row figwidth/5 50],...
    'Callback',{@era_viewprefs,'era_prefs',era_prefs,'era_data',era_data,...
    'inputs',inputs}); 

%tag gui
era_gui.Tag = 'era_gui';

end

function era_svb(varargin)
%back button. takes user back to era_start

%close gui
close(varargin{3});

%go back to era_start
era_start;

end

function era_svh(varargin)
%parses inputs to era_relfigures for displaying figures
%
%Input
% era_data - ERA Toolbox data structure array
% era_prefs - ERA Toolbox preferences structure array
%

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inputs
ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

%check whether the dependability estimate provided is numeric and between 0
%and 1
depeval = depcheck(str2double(inputs.h(1).String));

era_prefs.view.depvalue = str2double(inputs.h(1).String);
era_prefs.view.plotdep = inputs.h(2).Value;
era_prefs.view.ploticc = inputs.h(3).Value;
era_prefs.view.inctrltable = inputs.h(4).Value;
era_prefs.view.overalltable = inputs.h(5).Value;
era_prefs.view.showstddevt = inputs.h(6).Value;
era_prefs.view.showstddevf = inputs.h(7).Value;

%if the dependability estimate was not numeric or between 0 and 1, give the
%user an error and take the user back.
if depeval ~= 0 
    
    %check if era_gui is open.
    era_gui = findobj('Tag','era_gui');
    if ~isempty(era_gui)
        pos = era_gui.Position;
        close(era_gui);
    else
        pos=[400 400 550 550];
    end
    
    %create error text
    errorstr = {};
    errorstr{end+1} = 'The dependability estimate must be numeric';
    errorstr{end+1} = 'and between 0 and 1 (inclusive)';
    
    %display error prompt
    errordlg(errorstr);
    
    %execute era_startview_fig with the new preferences
    era_startview_fig('era_prefs',era_prefs,'era_data',era_data);
    
    return;
end

%pass inputs from gui to era_relfigures
era_relfigures('era_data',era_data,...
    'depcutoff',era_prefs.view.depvalue,...
    'plotdep',era_prefs.view.plotdep,...
    'ploticc',era_prefs.view.ploticc,...
    'showinct',era_prefs.view.inctrltable,...
    'showoverallt',era_prefs.view.overalltable,...
    'showstddevt',era_prefs.view.showstddevt,...
    'plotbetstddev',era_prefs.view.showstddevf,...
    'plotdepline',era_prefs.view.plotdepline,...
    'plotntrials',era_prefs.view.ntrials,...
    'meascutoff',era_prefs.view.meascutoff,...
    'depcentmeas',era_prefs.view.depcentmeas);

end


function era_viewprefs(varargin)
%displays various preferences for plotting or summarizing data
%
%Input
% era_data - ERA Toolbox data structure array
% era_prefs - ERA Toolbox preferences structure array

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%find inputs
ind = find(strcmp('inputs',varargin),1);
inputs = varargin{ind+1};

%check whether the dependability estimate provided is numeric and between 0
%and 1
depeval = depcheck(str2double(inputs.h(1).String));

%if the dependability estimate was not numeric or between 0 and 1, give the
%user an error and take the user back.
if depeval ~= 0 
    %create error text
    errorstr = {};
    errorstr{end+1} = 'The dependability estimate must be numeric';
    errorstr{end+1} = 'and between 0 and 1 (inclusive)';
    
    %display error prompt
    errordlg(errorstr);
    
    %execute era_startview_fig with the new preferences
    era_startview_fig(h_view_gui.filename,h_view_gui.pathname,'inputs',...
        h_view_gui.inputs,'viewprefs',initialprefs);
    
    return;
end

era_prefs.view.depvalue = str2double(inputs.h(1).String);
era_prefs.view.plotdep = inputs.h(2).Value;
era_prefs.view.ploticc = inputs.h(3).Value;
era_prefs.view.inctrltable = inputs.h(4).Value;
era_prefs.view.overalltable = inputs.h(5).Value;
era_prefs.view.showstddevt = inputs.h(6).Value;
era_prefs.view.showstddevf = inputs.h(7).Value;

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    pos = era_gui.Position;
    close(era_gui);
else
    pos=[400 400 550 550];
end

%define list for plotting dependability against number of trials
deplist = {'Lower Limit' 'Point Estimate' 'Upper Limit'};

%define list for central tendency measures
centlist = {'Mean' 'Median'};

%define space between rows and first row location
rowspace = 35;
row = pos(4) - rowspace*2;

%define locations of column 1 and 2 for the gui
lcol = 30;
rcol = (pos(3)/2+20);

%create the basic era_prefs
era_gui = figure('unit','pix',...
  'position',pos,...
  'menub','no',...
  'name','Specify Processing Preferences',...
  'numbertitle','off',...
  'resize','off');    

%print the gui headers
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [pos(4)/8 row pos(4)/3 25]);  

uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[4.4*pos(4)/8 row pos(4)/3 25]);

%next row
row = row - rowspace*2;

%which lines should be plotted on depplot
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Lines to plot for dependability',...
    'Position', [lcol row pos(4)/2 35]);  

newprefs.plotdepline = uicontrol(era_gui,'Style','listbox',...
    'fontsize',era_prefs.guis.fsize,...
    'String',deplist,'Min',1,'Max',1,'Value',era_prefs.view.plotdepline,...
    'Position', [rcol row pos(4)/3 50]);  

%next row
row = row - rowspace*2;

%number of trials to plot
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Number of trials to plot for dependability estimates',...
    'Position', [lcol row+5 pos(4)/2 35]);  

newprefs.ntrials = uicontrol(era_gui,...
    'Style','edit','fontsize',era_prefs.guis.fsize,...
    'String',era_prefs.view.ntrials,... 
    'Position', [rcol row+21 pos(4)/3 25]);  

%next row
row = row - rowspace*2.2;

%how to determine cutoff
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Estimate to use for trial cutoffs',...
    'Position', [lcol row pos(4)/2 35]);  

newprefs.meascutoff = uicontrol(era_gui,...
    'Style','listbox','fontsize',era_prefs.guis.fsize,...
    'String',deplist,'Min',1,'Max',1,...
    'Value',era_prefs.view.meascutoff,...
    'Position', [rcol row pos(4)/3 50]);  

%next row
row = row - rowspace*2.2;

%measure of central tendendcy for overall dependability
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    'Measure of central tendency for overall dependability calculations',...
    'Position', [lcol row pos(4)/2 35]);  

newprefs.depcentmeas = uicontrol(era_gui,'Style','listbox',...
    'fontsize',era_prefs.guis.fsize,...
    'String',centlist,'Min',1,'Max',1,...
    'Value',era_prefs.view.depcentmeas,...
    'Position', [rcol row pos(4)/3 40]);  

%next row with extra space
row = row - rowspace*2.5;

%Create a back button that will save inputs for preferences
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Save',...
    'Position', [pos(4)/8 row pos(4)/3 40],...
    'Callback',{@era_prefs_save,'era_prefs',era_prefs,'era_data',...
    era_data,'newprefs',newprefs}); 

%Create button that will go back to era_gui without saving
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [4.4*pos(4)/8 row pos(4)/3 40],...
    'Callback',{@era_prefs_back,'era_prefs',era_prefs,'era_data',...
    era_data});

%tag gui
era_gui.Tag = 'era_gui';

end

function era_prefs_back(varargin)
%if the back button was pressed the inputs will not be saved

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

%execute era_startview_fig with the old preferences
era_startview_fig('era_prefs',era_prefs,'era_data',era_data);

end

function era_prefs_save(varargin)
%if the save button was pressed use new inputs

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%find newprefs
ind = find(strcmp('newprefs',varargin),1);
newprefs = varargin{ind+1};

%pull new preferences
era_prefs.view.plotdepline = newprefs.plotdepline.Value;
era_prefs.view.ntrials = str2double(newprefs.ntrials.String);
era_prefs.view.meascutoff = newprefs.meascutoff.Value;
era_prefs.view.depcentmeas = newprefs.depcentmeas.Value;

%check if era_gui is open
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

%execute era_startview_fig with the new preferences
era_startview_fig('era_prefs',era_prefs,'era_data',era_data);

end

function checkout = depcheck(depvalue)
%ensure that the provided dependability estimate is numeric and between 0
%and 1
%
%Input
% depvalue - dependability threshold estimate from era_startview_fig
%
%Output
% checkout
%   0: dependability estimate is numeric and between 0 and 1
%   1: dependability is string
%   2: dependability is not between 0 and 1

%check whether depvalue is numeric
if isnan(depvalue)
    checkout = 1;
else
    %check whether depvalue is between 0 and 1
    if depvalue > 0 && depvalue < 1
        checkout = 0;
    else
        checkout = 2;
    end
end

end