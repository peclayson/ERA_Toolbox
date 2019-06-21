function era_startview_trt(varargin)
%Prepares the data from multiple occasions for viewing and lets user 
% specify which tables and figures to present
%
%era_startview_trt('era_prefs',era_prefs,'era_data',era_data)
%
%Last Updated 6/21/19
%
%Required Input
% era_data - ERA Toolbox data structure array
% era_prefs - ERA Toolbox preferences structure array
%
%Output
% No variables will be outputted to the Matlab workspace. Based on the
%  inputs from this gui, era_relfigures will be executed to display various
%  figures and tables (for more information about the tables and figures
%  see the user manual for the ERA toolbox)

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
% by Peter Clayson (6/21/19)
% peter.clayson@gmail.com
%

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

%define parameters for figure position
figwidth = 550;
figheight = 650;

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

%check that a default fsize has been defined
if ~isfield(era_prefs,'fsize')
    era_prefs.guis.fsize = get(0,'DefaultTextFontSize');
end

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
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

%next row
row = row - (rowspace*.45);

%Print the name of the measurement analyzed
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Measurement:  ' era_data.proc.measheader],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]); 

%next row
row = row - (rowspace*1.4);

str = sprintf(['The reliability threshold to use for retaining data\n'...
    'Participants that do not have enough trials to meet this reliability\n'...
    'threshold will be recommended for exclusion']);

%Print the text for reliability cutoff with a box for the user to specify
%the input
uicontrol(era_gui,...
    'Style','text',...
    'fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Reliability Cutoff:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.relcutoff = uicontrol(era_gui,...
    'Style','edit',...
    'fontsize',era_prefs.guis.fsize,...
    'String',era_prefs.view.relvalue,...
    'Position', [rcol+5 row+6 figwidth/4 25]);  

%next row
row = row - rowspace-5;

str = sprintf(['Choose dependability to use the absolute error variance\n'...
    'Choose generalizability to use the relative error variance\n'...
    'See user manual for more information about coefficients']);

%Provide the user with the option to choose between dependability and
%generalizability coefficients
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Type of G-Theory Coefficient:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.depgen = uicontrol(era_gui,...
    'Style','pop',...
    'fontsize',era_prefs.guis.fsize,...
    'String',{'Dependability' 'Generalizability'},...
    'Value',era_prefs.view.gcoeff,...
    'Position', [rcol row figwidth/4 25]);  

%next row
row = row - rowspace-10;

str = sprintf(['For internal consistency, use a coefficient of equivalence\n'...
    'For test-retest reliability, use a coefficient of stability\n'...
    'See user manual for more information about coefficients']);

%Provide the user with the option to choose between coefficients of 
%equivalence and coefficients of stability
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Type of Reliability Coefficient:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);  

inputs.equistab = uicontrol(era_gui,...
    'Style','pop',...
    'fontsize',era_prefs.guis.fsize,...
    'String',{'Equivalence' 'Stability'},...
    'Value',era_prefs.view.reltype,...
    'Position', [rcol row figwidth/4 25]);  

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

chckstr = 'Checked = Yes; Unchecked = No';
str = sprintf(['Display a plot that shows the impact of the number of\n'...
    'trials retained for averaging on reliability estimates']);

%reliability with increasing trials
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: Number of Trials v Reliability',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.plotrel = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.plotrel,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%plot ICCs
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: Intraclass Correlation Coefficients',...
    'Tooltip','Display a plot of the intraclass correlation coefficients',...
    'Position', [lcol row figwidth/2 40]);  

inputs.ploticc = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.ploticc,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%plot between-person standard deviations
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Plot: Between-Person Standard Deviations',...
    'Tooltip','Display a plot showing the between-person standard deviations',...
    'Position', [lcol row figwidth/2 40]);  

inputs.plotbetsd = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.showstddevf,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row
row = row - rowspace;

str = sprintf(['Display a table showing the number of trials needed\n',...
    'to obtain the specified reliability threshold and the\n',...
    'reliability point estimate and credible interval for the trial cutoff']);
    
%reliability cutoff table
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Trial Cutoffs for Specified Reliability Threshold',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.relcutt = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.inctrltable,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

str = sprintf(['Display a table that summarizes the number of participants\n'...
    'with data that satisfy the reliability threshold, the number of\n'...
    'participants without data the satisfy the reliability threshold\n',...
    'the overall reliability point estimate and credible interval,\n'...
    'and the trial summary information (min, mean, median, max)']);


%overall reliability table
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Overall Reliability and Summary Information',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.overallt = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.overalltable,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]); 

%next row
row = row - rowspace;

%between- and within-person standard deviation tables
str = sprintf(['Display a table that shows the between- and within-person\n',...
    'standard deviations and intraclass correlation coefficients']);

uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Sources of Variance',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.sdt = uicontrol(era_gui,'Style','checkbox',...
    'Value',era_prefs.view.showstddevt,...
    'Tooltip',chckstr,...
    'Position', [rcol row+20 figwidth/2 25]);

%next row for buttons
row = row - rowspace*1.4;

%Create a back button that will take the user back to era_start
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [.25*figwidth/9 row figwidth/6 50],...
    'Callback',{@era_svb,era_gui});

%Create a button to let the user load a different .erat file
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Load New File',...
    'Position', [2*figwidth/9 row figwidth/6 50],...
    'Callback',{@era_view_loadnewfile,era_gui});

%Create button that will check the inputs and begin processing the data
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [3.75*figwidth/9 row figwidth/6 50],...
    'Callback',{@era_svh,'era_prefs',era_prefs,'era_data',era_data,...
    'inputs',inputs}); 

%Create button that will display preferences
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Preferences',...
    'Position', [5.5*figwidth/9 row figwidth/6 50],...
    'Callback',{@era_viewprefs,'era_prefs',era_prefs,'era_data',era_data,...
    'inputs',inputs}); 

%Create button that will close any open figures other than this gui
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Close Windows',...
    'Position', [7.25*figwidth/9 row figwidth/6 50],...
    'Tooltip','Close open Matlab windows other than this one',...
    'Callback',{@era_closefigs}); 

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

function era_view_loadnewfile(varargin)
%give user the chance to load a new file

%close gui
close(varargin{3});
era_startview;

end

function era_closefigs(varargin)
%button to close Matlab windows other than the era_gui

%get the hanles for all open figures
currfigs = findall(0,'Type','figure');

%close all of the figures that are not the Specify Inputs gui
for i = 1:length(currfigs)
    if ~strcmp(currfigs(i).Name,'Specify Inputs')
        close(currfigs(i));
    end
end


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

%check whether the reliability estimate provided is numeric and between 0
%and 1
releval = relcheck(str2double(inputs.relcutoff.String));

era_prefs.view.relvalue = str2double(inputs.relcutoff.String);
era_prefs.view.gcoeff = inputs.depgen.Value;
era_prefs.view.reltype = inputs.equistab.Value;
era_prefs.view.plotrel = inputs.plotrel.Value;
era_prefs.view.ploticc = inputs.ploticc.Value;
era_prefs.view.inctrltable = inputs.relcutt.Value;
era_prefs.view.overalltable = inputs.overallt.Value;
era_prefs.view.showstddevt = inputs.sdt.Value;
era_prefs.view.showstddevf = inputs.plotbetsd.Value;

%if the reliability estimate was not numeric or between 0 and 1, give the
%user an error and take the user back.
if releval ~= 0 
    
    %check if era_gui is open.
    era_gui = findobj('Tag','era_gui');
    if ~isempty(era_gui)
        close(era_gui);
    end
    
    %create error text
    errorstr = {};
    errorstr{end+1} = 'The reliability estimate must be numeric';
    errorstr{end+1} = 'and between 0 and 1 (inclusive)';
    
    %display error prompt
    errordlg(errorstr);
    
    %execute era_startview_sing with the new preferences
    era_startview_trt('era_prefs',era_prefs,'era_data',era_data);
    
    return;
end

%pass inputs from gui to era_relfigures
era_relfigures('era_data',era_data,...
    'relcutoff',era_prefs.view.relvalue,...
    'plotrel',era_prefs.view.plotrel,...
    'ploticc',era_prefs.view.ploticc,...
    'showinct',era_prefs.view.inctrltable,...
    'showoverallt',era_prefs.view.overalltable,...
    'showstddevt',era_prefs.view.showstddevt,...
    'plotbetstddev',era_prefs.view.showstddevf,...
    'plotrelline',era_prefs.view.plotrelline,...
    'plotntrials',era_prefs.view.ntrials,...
    'meascutoff',era_prefs.view.meascutoff,...
    'relcentmeas',era_prefs.view.relcentmeas);

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

if ~isempty(ind)
    inputs = varargin{ind+1};

    %check whether the reliability estimate provided is numeric and between 0
    %and 1
    releval = relcheck(str2double(inputs.h(1).String));

    %if the reliability estimate was not numeric or between 0 and 1, give the
    %user an error and take the user back.
    if releval ~= 0 
        %create error text
        errorstr = {};
        errorstr{end+1} = 'The reliability estimate must be numeric';
        errorstr{end+1} = 'and between 0 and 1 (inclusive)';

        %display error prompt
        errordlg(errorstr);

        %execute era_startview_sing with the new preferences
        era_startview_sing(h_view_gui.filename,h_view_gui.pathname,'inputs',...
            h_view_gui.inputs,'viewprefs',initialprefs);
        
        return;
    end
    
    era_prefs.view.relvalue = str2double(inputs.relcutoff.String);
    era_prefs.view.gcoeff = inputs.depgen.Value;
    era_prefs.view.reltype = inputs.equistab.Value;
    era_prefs.view.plotrel = inputs.plotrel.Value;
    era_prefs.view.ploticc = inputs.ploticc.Value;
    era_prefs.view.inctrltable = inputs.relcutt.Value;
    era_prefs.view.overalltable = inputs.overallt.Value;
    era_prefs.view.showstddevt = inputs.sdt.Value;
    era_prefs.view.showstddevf = inputs.plotbetsd.Value;
end

%check if era_gui is open.
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    pos = era_gui.Position;
    close(era_gui);
else
    pos=[400 400 550 650];
end

%define list for plotting reliability against number of trials
rellist = {'Lower Limit' 'Point Estimate' 'Upper Limit'};

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

str = sprintf(['Indicate the estimate that should be plotted in the\n',...
    'figure that shows the relationship between the number of trials\n',...
    'retained for averaging and reliability']);

%which lines should be plotted on relplot
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Line to plot for reliability',...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/2 35]);  

str = sprintf(['Lower Limit of Credible Interval\n',...
    'Reliability Point Estimate\n',...
    'Upper Limit of Credible Interval']);

newprefs.plotrelline = uicontrol(era_gui,'Style','listbox',...
    'fontsize',era_prefs.guis.fsize,...
    'String',rellist,'Min',1,'Max',1,'Value',era_prefs.view.plotrelline,...
    'Tooltip',str,... 
    'Position', [rcol row pos(4)/3 50]);  

%next row
row = row - rowspace*2;

str = sprintf(['Indicate the number of trials that should be plotted in the\n',...
    'figure that shows the relationship between the number of trials\n',...
    'retained for averaging and reliability']);

%number of trials to plot
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Number of trials to plot',...
    'Tooltip',str,...
    'Position', [lcol row+5 pos(4)/2 35]);  

newprefs.ntrials = uicontrol(era_gui,...
    'Style','edit','fontsize',era_prefs.guis.fsize,...
    'String',era_prefs.view.ntrials,... 
    'Position', [rcol row+21 pos(4)/3 25]);  

%next row
row = row - rowspace*2.2;

str = sprintf(['Indicate which reliability estimate should be used\n',...
    'for the reliability threshold that deems data as reliable']);

%how to determine cutoff
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Estimate to use for trial cutoffs',...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/2 35]);  

str = sprintf(['Lower Limit of Credible Interval\n',...
    'Reliability Point Estimate\n',...
    'Upper Limit of Credible Interval']);

newprefs.meascutoff = uicontrol(era_gui,...
    'Style','listbox','fontsize',era_prefs.guis.fsize,...
    'String',rellist,'Min',1,'Max',1,...
    'Value',era_prefs.view.meascutoff,...
    'Tooltip',str,...
    'Position', [rcol row pos(4)/3 50]);  

%next row
row = row - rowspace*2.2;

str = sprintf(['Indicate the measure of central tendency to use for the\n',...
    'estimation of the overall reliability of the dataset']);

%measure of central tendendcy for overall reliability
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',...
    'Measure of central tendency for overall reliability calculations',...
    'Tooltip',str,...
    'Position', [lcol row pos(4)/2 35]);  

newprefs.relcentmeas = uicontrol(era_gui,'Style','listbox',...
    'fontsize',era_prefs.guis.fsize,...
    'String',centlist,'Min',1,'Max',1,...
    'Value',era_prefs.view.relcentmeas,...
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

%execute era_startview_trt with the old preferences
era_startview_trt('era_prefs',era_prefs,'era_data',era_data);

end

function era_prefs_save(varargin)
%if the save button was pressed use new inputs

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%find newprefs
ind = find(strcmp('newprefs',varargin),1);
newprefs = varargin{ind+1};

%pull new preferences
era_prefs.view.plotrelline = newprefs.plotrelline.Value;
era_prefs.view.ntrials = str2double(newprefs.ntrials.String);
era_prefs.view.meascutoff = newprefs.meascutoff.Value;
era_prefs.view.relcentmeas = newprefs.relcentmeas.Value;

%check if era_gui is open
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end

if era_prefs.view.ntrials > 1
    %execute era_startview_trt with the new preferences
    era_startview_trt('era_prefs',era_prefs,'era_data',era_data);
end

%make sure the user has defined at least 2 trials to plot for the figure
if era_prefs.view.ntrials <= 1 
    errordlg(['Please specify at least two trials to plot ', ...
        'for reliability estimates']);
    era_prefs.view.ntrials = 50;

    era_viewprefs('era_prefs',era_prefs,'era_data',era_data);
    
end

end

function checkout = relcheck(relvalue)
%ensure that the provided reliability estimate is numeric and between 0
%and 1
%
%Input
% relvalue - reliability threshold estimate from era_startview_sing
%
%Output
% checkout
%   0: reliability estimate is numeric and between 0 and 1
%   1: reliability is string
%   2: reliability is not between 0 and 1

%check whether relvalue is numeric
if isnan(relvalue)
    checkout = 1;
else
    %check whether relvalue is between 0 and 1
    if relvalue > 0 && relvalue < 1
        checkout = 0;
    else
        checkout = 2;
    end
end

end