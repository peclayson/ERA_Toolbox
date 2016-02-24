function ra_start
%
%Initiate Matlab gui to use the ERP Reliability Analysis (RA) toolbox
%
%version .80 - Last Updated 2/22/16
%

%The RA toolbox uses generalizability theory as a method for evaluating 
%realibity of ERP data. Dependability estimates (generalizability-theory
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
% Pe. Psychophysiology, 52(6), 790?800. http://doi.org/10.1111/psyp.12401
%
%The notion of reporting estimates of reliability in all ERP studies and 
% this toolbox are specifically discussed in 
%
% <insert citation here>
%
%Please cite both papers when using the RA toolbox (Baldwin et al. paper
%for the formulas and concept; Clayson and Miller paper for the toolbox)
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
%Matlab files contained in RA toolbox
%
% Gui-related files
%  ra_startproc - gui to specify inputs for the analysis of data using 
%   cmdstand (runs ra_loadfile and ra_computerel)
%  ra_startview - gui to specify inputs for displaying processed data
%   (runs ra_relfigures)
%
% Data analysis files
%  ra_loadfile - function to load data files and prepare them for 
%   processing
%  ra_computerel - function to run data through cmdstan for dependability
%   analyses
%  ra_relfigures - function to display information about dependability
%   (figures and tables)
%
%Required dependents not contained in download of RA toolbox
% MatlabStan and its dependents: CmdStan and MatlabProcessManager
% Instructions can be found at 
%  https://github.com/brian-lau/MatlabStan/wiki/Getting-Started
% These dependents are necessary for computing dependability data. IF you
%  are only interested in viewing already processed data from the RA 
%  toolbox then MatlabStan, CmdStan, and MatlabProcessManager are not
%  necessary.

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
% by Peter Clayson (2/22/16)
% peter.clayson@gmail.com
%

%Output info about RA Toolbox
fprintf('\n\n\nERP Reliability Analysis Toolbox Version .80\n\n');


%check whether dependencies are contained in the Matlab path
%first look for RA toolbox files 
if exist('ra_startproc.m','file') ~= 2 || ...
        exist('ra_startview.m','file') ~= 2 || ...
        exist('ra_computerel.m','file') ~= 2 || ...
        exist('ra_loadfile.m','file') ~= 2 || ...
        exist('ra_relfigures.m','file') ~= 2 

    dlg = {'Warning: Dependencies for the RA toolbox are not located';...
        'in the Matlab path. Please include the folder containing the';...
        'scripts for the RA toolbox in your Matlab path.';...
        'Please see RA_UserManual.pdf for more information'};
    
    for i = 1:length(dlg)
        fprintf('%s\n',dlg{i});
    end
    fprintf('\n\n');
    return;
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
    

    dlg = {'Warning: Dependencies for the RA toolbox are not located';...
        'in the Matlab path. Please include the folders containing the';...
        'scripts for CmdStan, MatlabProcessManager, and MatlabStan';...
        'in your Matlab path.';...
        'Please see RA_UserManual.pdf for more information'};
    
    for i = 1:length(dlg)
        fprintf('%s\n',dlg{i});
    end
    fprintf('\n\n');
    return;
end

%set the gui for indicating whether the user wants to process data or view
%previously processed data

%define parameters for figure position
figwidth = 400;
figheight = 200;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%initialize gui
ra_gui= figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');
movegui(ra_gui,'center');

%Print the name of the loaded dataset
uicontrol(ra_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','What would you like to do?',...
    'Position',[0 row figwidth 25]);          

%Create a button that will take the user to the gui for setting the inputs
%to process data
uicontrol(ra_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','<html><center>Process <br>New Data',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',@ra_startproc); 

%Create button that will take the user to the gui for setting the inputs
%for viewing the data
uicontrol(ra_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','View Results',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@ra_startview}); 

%display gui
set(ra_gui,'Visible','on');

end
