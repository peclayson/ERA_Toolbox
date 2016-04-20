function era_start
%
%Initiate Matlab gui to use the ERP Reliability Analysis (ERA) toolbox
%
%version 0.3.1 - Last Updated 4/18/16
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
% Baldwin, S. A., Larson, M. J., & Clayson, P. E. (2015). The dependability 
% of electrophysiological measurements of performance monitoring in a 
% clinical sample: A generalizability and decision analysis of the ERN and 
% Pe. Psychophysiology, 52(6), 790-800. http://doi.org/10.1111/psyp.12401
%
%The notion of reporting estimates of reliability in all ERP studies and 
% this toolbox are specifically discussed in 
%
% Clayson, P. E., & Miller, G. A. (under review). Psycometric
% considerations in the measurement of event-related brain potentials:
% Guidelines for measurement and reporting
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
%Matlab files contained in ERA toolbox
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
%
%Required dependents not contained in download of ERA toolbox
% MatlabStan and its dependents: CmdStan and MatlabProcessManager
% Instructions can be found at 
%  https://github.com/brian-lau/MatlabStan/wiki/Getting-Started
% These dependents are necessary for computing dependability data. 
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
% added subroutines to dependents check

%set version number of ERA Toolbox
eraver = '0.3.1';

%Output info about ERA Toolbox
fprintf('\n\n\nERP Reliability Analysis Toolbox Version %s\n\n',eraver);

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
        exist('era_reruncheck.m','file') ~= 2

    dlg = {'Warning: Dependencies for the ERA toolbox are not located';...
        'in the Matlab path. Please include the folder containing the';...
        'scripts for the RA toolbox in your Matlab path.';...
        'Please see RA_UserManual.pdf for more information'};
    
    for i = 1:length(dlg)
        fprintf('%s\n',dlg{i});
    end
    fprintf('\n\n');
    return;
    
else
    fprintf('\nERA Toolbox files found');
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
        'Please see RA_UserManual.pdf for more information'};
    
    for i = 1:length(dlg)
        fprintf('%s\n',dlg{i});
    end
    fprintf('\n\n');
    return;
else
    fprintf('\nCmdStan, MatlabProcessManager, and MatlabStan files found');
end

fprintf('\n\n');

%set the gui for indicating whether the user wants to process data or view
%previously processed data

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

%Write text
uicontrol(era_gui,'Style','text','fontsize',fsize+2,...
    'HorizontalAlignment','center',...
    'String','What would you like to do?',...
    'Position',[0 row figwidth 25]);          

%Create a button that will take the user to the gui for setting the inputs
%to process data
uicontrol(era_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','<html><center>Process <br>New Data',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_startproc,era_gui}); 

%Create button that will take the user to the gui for setting the inputs
%for viewing the data
uicontrol(era_gui,'Style','push','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','View Results',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_startview,era_gui}); 

%tag gui
era_gui.Tag = 'era_gui';

%display gui
set(era_gui,'Visible','on');

end
