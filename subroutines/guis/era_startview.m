function era_startview(varargin)
%Loads data for viewing and determines whether to load gui for examining
% single occasion or multiple occasion data
%
%era_startview('file','/Users/REL/SomeERAData.mat')
%
%Last Updated 6/21/19
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
%
%7/26/16 PC
% added check to make sure that at least 2 trials were requested for the
%  number of trials and dependability plot
%
%7/27/16 PC
% got rid of some code that was no longer used
%
%9/18/16 PC
% added tooltips (text that appears when you hover over a gui
%  property)
% changed the text that is displayed in the gui to be more concise, since
%  additional explanation is provided in tooltip
% added a button to close all open figures other than Specify Inputs gui
%
%1/19/17 PC
% updated copyright
%
%8/16/17 PC
% fixed bug with era_startview not working correctly unless 
%  era_prefs.guis.fsize had already been defined
%
%8/23/17 PC
% added a button for loading a new file from the gui
%
%6/22/18 PC
% added which measurement was processed
%
%6/17/19 PC
% fixed bug: era_data not loaded when specified in era_startview
%
%6/21/19 PC
% this function was turned into a wrapper to determine whether data
%  contain just one occasion (just internal consistency) or mulitple
%  occasions (internal consistency + trt)

%somersault through varargin inputs to check for era_prefs and era_data
[era_prefs, era_data] = era_findprefsdata(varargin);

%see if the file for the figures and tables has been specified in
%varargin
if ~isempty(varargin) && (isempty(era_data) && isempty(era_prefs))
    
    %check if data file has been provided
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = varargin{ind+1};
        
        %if file has been provided, load it
        load(file,'-mat');
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
%pic between two possible guis: gui for single session data and a gui for
%multiple session data
if ~isfield(era_data.rel,'analysis') ||... 
        (isfield(era_data.rel,'analysis') && ...
        strcmp(era_data.rel.analysis,'ic'))
    era_startview_sing('era_prefs',era_prefs,'era_data',era_data);
elseif isfield(era_data.rel,'analysis') &&... 
        strcmp(era_data.rel.analysis,'trt')
    era_startview_trt('era_prefs',era_prefs,'era_data',era_data);
end

end
