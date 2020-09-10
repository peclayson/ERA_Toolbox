function era_start
%
%Initiate Matlab gui to use the ERP Reliability Analysis (ERA) toolbox
%
%version 0.5.1 - Last Updated 8/31/20
%

%The ERA toolbox uses generalizability theory as a method for evaluating 
%reliability of ERP data. Dependability estimates (generalizability-theory
%analog to reliability) can be computed for any number of groups or events.
%The influence of number of trials on dependability of measurements will
%also be determined, and a recommended cutoff for inclusion of ERP data 
%will be provided based on the stability of measurement as the number of 
%trials included in a single-subject average for a given group and event 
%increases.
%
%A description of how to apply generalizability theory to ERP data can be
% found in
%
% Clayson, P. E., & Miller, G. A. (2017). ERP Reliability Analysis
% (ERA) Toolbox: An open-source toolbox for analyzing the reliability of
% event-related potentials. International Journal of Psychophysiology, 111,
% 68-79. doi: 10.1016/j.ijpsycho.2016.10.012
%
% Baldwin, S. A., Larson, M. J., & Clayson, P. E. (2015). The dependability 
% of electrophysiological measurements of performance monitoring in a 
% clinical sample: A generalizability and decision analysis of the ERN and 
% Pe. Psychophysiology, 52, 790-800. doi: 10.1111/psyp.12401
%
% Clayson, P. E., Carbine, K. A., Baldwin, S. A., Olsen, J. A., & 
% Larson, M. J. (under review). Using generalizability theory and the ERP 
% Reliability Analysis (ERA) Toolbox for assessing test-retest reliability 
% of ERP scores: Algorithms, framework, and implementation. 
%
% Clayson, P. E., Brush, C. J., & Hajcak, G. (under review). Data quality
% and reliability metrics for event-related potentials (ERPs): The utility 
% of subject-level reliability.
%
%The notion of reporting estimates of reliability in all ERP studies and 
% this toolbox are specifically discussed in 
%
% Clayson, P. E., & Miller, G. A. (2017). Psycometric
% considerations in the measurement of event-related brain potentials:
% Guidelines for measurement and reporting. International Journal of
% Psychophysiology, 111, 57-67. doi: 10.1016/j.ijpsycho.2016.09.005
%
%
%Input
% There are no required inputs to execute this script. 
% However, paths to the data to be processed or viewed will need to be 
%  provided.
%
%Output
% This script will not output any variables to the workspace. However, if
%  the user chooses, various files (data, tables, plots) could be saved
%  from the guis initiated by this script.
%
%High-level matlab files contained in ERA toolbox
%
% Gui-related files
%  era_startproc - gui to specify inputs for the analysis of data using 
%   cmdstand (runs era_loadfile and era_computerel)
%  era_startview - gui to specify inputs for displaying processed data
%   (runs era_relfigures)
%
% Data analysis files
%  era_loadfile - function to load data files and prepare them for 
%   processing
%  era_computerel - function to run data through cmdstan for dependability
%   analyses
%  era_relfigures - function to display information about dependability
%   (figures and tables)
%  era_defaults - function to define default settings for processing and
%   viewing data
%
%Required dependents not contained in download of ERA toolbox
% MatlabStan and its dependents: CmdStan and MatlabProcessManager
% Instructions can be found at 
%  https://github.com/brian-lau/MatlabStan/wiki/Getting-Started
% These dependents are necessary for computing dependability data. 
%

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

%If you have any questions or comments of if you run into any bugs, please
%contact me at peter.clayson@gmail.com. I will do my best to respond
%promptly. Please note that if you run into a bug/error I have likely
%not encountered it, but I will do my best to provide any support
%and help I can to fix it.

%History 
% by Peter Clayson (4/18/16)
% peter.clayson@gmail.com
%
%4/20/16 PC
% added subroutine to check whether there is a new release posted on github
%
%7/21/16 PC
% changes associated with adding a new data structure
% change with setting version number
%
%8/6/16 PC
% ERA Toolbox directories automatically added to path
% automatic installation function added
% added gui to ask whether the dependents should be automatically installed
%
%8/21/16 PC
% added check to make sure subroutines/plotting is on path
%
%9/18/16 PC
% added return when checking for installation so function does not continue
%
%1/19/17 PC
% updated citations, version number, and copyright
%
%6/23/17 PC
% updated and tested new versions of dependents
%
%6/24/17 PC
% added functions to update version of ERA_Toolbox to current version
% added functions to update old versions of dependents (except for
%  MatlabStan)
%
%6/25/17 PC
% change dimensions of gui that asks whether the dependents should be
%  updated
%
%9/7/17 PC
% fixed placement of text in era_ask2updatedependents

%check whether dependencies are contained in the Matlab path
%first look for ERA toolbox files
fprintf('\nEnsuring dependents are found in the Matlab path\n');
if exist('era_startproc.m','file') ~= 2 || ...
        exist('era_startview.m','file') ~= 2 || ...
        exist('era_computerel.m','file') ~= 2 || ...
        exist('era_loadfile.m','file') ~= 2 || ...
        exist('era_relfigures.m','file') ~= 2 || ...
        exist('era_checkconv.m','file') ~= 2 || ...
        exist('era_readtable.m','file') ~= 2 || ...
        exist('era_reruncheck.m','file') ~= 2 || ...
        exist('era_updatecheck.m','file') ~= 2 || ...
        exist('era_defaults.m','file') ~= 2 || ...
        exist('era_depvtrialsplot.m','file') ~= 2 || ...
        exist('era_ptintervalplot.m','file') ~= 2 || ...
        exist('era_installdependents.m','file') ~= 2 || ...
        exist('era_checkversionsofdeps.m','file') ~= 2
    
    %find where the era files are located and add the directory and
    %sub-directories
    eras_path = which('era_start');
    udir = fileparts(eras_path);
    %add files to the Matlab path and save it
    addpath(genpath(udir));
    savepath;
    fprintf('Directories for ERA Toolbox files added to Matlab path\n');
    
else
    fprintf('ERA Toolbox files found\n');
end

%look for Stan dependents (the important files, assuming the rest are in
%the Matlab path
if exist('mcmc.m','file') ~= 2 || ...
        exist('stan.m','file') ~= 2 || ...
        exist('StanFit.m','file') ~= 2 || ...
        exist('StanModel.m','file') ~= 2 || ...
        exist('processManager.m','file') ~= 2 || ...
        exist('processState.m','file') ~= 2 || ...
        exist('makefile','file') ~= 2 || ...
        exist('test-all.sh','file') ~= 2 || ...
        exist('runCmdStanTests.py','file') ~= 2 
    

    dlg = {'Warning: Dependencies for the ERA toolbox are not located';...
        'in the Matlab path. Please include the folders containing the';...
        'scripts for CmdStan, MatlabProcessManager, and MatlabStan';...
        'in your Matlab path.';...
        'Please see UserManual.pdf for more information'};
    
    for i = 1:length(dlg)
        fprintf('%s\n',dlg{i});
    end
    fprintf('\n\n');
    era_ask2install;
    return;
    
else
    
%     fprintf('CmdStan, MatlabProcessManager, and MatlabStan files found\n');
%     
%     %check if the most up-to-date releases are being used for each of the
%     %dependents
%     
%     depvercheck = era_checkversionsofdeps;
%     
%     switch depvercheck.cmdstan
%         case 0
%             fprintf('There is a new version of CmdStan available\n');
%         case 1
%             fprintf('CmdStan version is current\n');
%         case 2
%             warning(strcat('Installed CmdStan version is newer than the',...
%                 ' version tested for the Toolbox.',...
%                 ' Toolbox may not work properly as a result'));
%     end
%     
%     switch depvercheck.matlabprocessmanager
%         case 0
%             fprintf('There is a new version of MatlabProcessManager available\n');
%         case 1
%             fprintf('MatlabProcessManager version is current\n');
%         case 2
%             warning(strcat('Installed MatlabProcessManager version is newer than the',...
%                 ' version tested for the Toolbox.',...
%                 ' Toolbox may not work properly as a result'));
%     end
%     
%     fprintf('Checking the installed version of MatlabStan is not currently supported\n');
    
end

%pull version for the ERA Toolbox
eraver = era_defineversion;

%Output info about ERA Toolbox
fprintf('\n\n\nERP Reliability Analysis Toolbox Version %s\n\n',eraver);

%check whether running the newest release of the toolbox
% eravercheck = era_updatecheck(eraver);

fprintf('\n');

%load default preferences for processing and viewing data
era_prefs = era_defaults;

%attach the current version number to era_prefs
era_prefs.ver = eraver;

%create the gui for indicating whether the user wants to process data or view
%previously processed data

%define parameters for figure position
figwidth = 400;
figheight = 200;
era_prefs.guis.fsize = get(0,'DefaultTextFontSize');

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
era_gui= figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(era_gui,'center');

%Write text
uicontrol(era_gui,'Style','text','fontsize',era_prefs.guis.fsize+2,...
    'HorizontalAlignment','center',...
    'String','What would you like to do?',...
    'Position',[0 row figwidth 25]);          

%Create a button that will take the user to the gui for setting the inputs
%to process data
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','<html><center>Process <br>New Data',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_startproc,era_gui,'era_prefs',era_prefs}); 

%Create button that will take the user to the gui for setting the inputs
%for viewing the data
uicontrol(era_gui,'Style','push','fontsize',era_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','View Results',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_startview,era_gui,'era_prefs',era_prefs}); 

%tag gui
era_gui.Tag = 'era_gui';

%display gui
set(era_gui,'Visible','on');

%if applicable, ask the user whether ERA Toolbox should be updated
if eravercheck == 0
    
    era_ask2updatetoolbox;

%if applicable, ask the user whether dependents should be updated
elseif depvercheck.cmdstan == 0 || depvercheck.matlabprocessmanager == 0   
    
    era_ask2updatedeps(depvercheck);
    
end

end

function era_ask2install
%gui to ask the user whether the ERA Toolbox dependents should be installed

%define parameters for figure position
figwidth = 400;
figheight = 200;
fsize = get(0,'DefaultTextFontSize');

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
era_gui= figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(era_gui,'center');

str = ['ERA Toolbox dependents were not located. Would you like'...
    ' to install the dependents?'];

%Write text
uicontrol(era_gui,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 40]);          

%Create a button that will take install the dependents
uicontrol(era_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','<html><center>Install <br>Dependents',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@installdeps}); 

%Create button that quit
uicontrol(era_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Exit',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@giveup}); 

%tag gui
era_gui.Tag = 'era_gui';

%display gui
set(era_gui,'Visible','on');

end

function era_ask2updatetoolbox
%gui to ask the user whether the ERA Toolbox dependents should be installed

%define parameters for figure position
figwidth = 400;
figheight = 200;
fsize = get(0,'DefaultTextFontSize');

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
era_gui_update = figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(era_gui_update,'center');

str = ['You are using an old version of the ERA Toolbox.'...
    ' It is recommended that you update the toolbox.'...
    ' Would you like to do so now?'];

%Write text
uicontrol(era_gui_update,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 40]);          

%Create a button that will take install the dependents
uicontrol(era_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Yes',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@updateera}); 

%Create button that quit
uicontrol(era_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','No',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@donotupdateera}); 

%tag gui
era_gui_update.Tag = 'era_gui_update';

%display gui
set(era_gui_update,'Visible','on');

end

function era_ask2updatedeps(depvercheck)
%gui to ask the user whether the ERA Toolbox dependents should be installed

%define parameters for figure position
figwidth = 500;
figheight = 200;
fsize = get(0,'DefaultTextFontSize');

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*4;

%initialize gui
era_gui_update = figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(era_gui_update,'center');

deps2update = '';

if depvercheck.cmdstan == 0
    deps2update = [deps2update ' CmdStan'];
elseif depvercheck.matlabstan == 0
    deps2update = [deps2update ' MatlabStan'];
elseif depvercheck.matlabprocessmanager == 0
    deps2update = [deps2update ' MatlabProcessManager'];
end

str = ['You are using old version(s) of' deps2update...
    '. It is recommended that you update the toolbox.'...
    ' Would you like to do so now? WARNING: Doing so will delete the '...
    'old directories for the dependents to avoid confusion.'];

%Write text
uicontrol(era_gui_update,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 75]);          

%Create a button that will take install the dependents
uicontrol(era_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Yes',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@updatedeps,depvercheck}); 

%Create button that quit
uicontrol(era_gui_update,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','No',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@donotupdatedeps}); 

%tag gui
era_gui_update.Tag = 'era_gui_update';

%display gui
set(era_gui_update,'Visible','on');

end

function donotupdateera(varargin)
%if the user does not what the ERA toolbox updated, just continue on
%check if era_gui_update is open.
era_gui_update = findobj('Tag','era_gui_update');

if ~isempty(era_gui_update)
    close(era_gui_update);
end

fprintf('\nERP Relaibility Analysis Toolbox will not be updated\n\n'); 

end

function donotupdatedeps(varargin)
%if the user does not what the ERA toolbox updated, just continue on
%check if era_gui_update is open.
era_gui_update = findobj('Tag','era_gui_update');

if ~isempty(era_gui_update)
    close(era_gui_update);
end

fprintf('\nDependents will not be updated\n\n'); 

end

function updateera(varargin)
%if the user wants the toolbox updated, begin updating
close all;
era_updateera;
end

function updatedeps(varargin)
%if the user wants the dependents updated, begin updating
close all;
depvercheck = varargin{3};
era_installdependents('depvercheck',depvercheck);
end

function giveup(varargin)
%if the user does not want the dependents installed, then close everything
close all;
end

function installdeps(varargin)
%if the user wants the dependents installed, begin the installation
close all;
era_installdependents;
end
