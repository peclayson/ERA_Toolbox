function era_startview(varargin)
%Prepares the data for viewing and lets user specify what tables and
%figures to present
%
%era_startview('file','/Users/REL/SomeERAData.mat')
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
% by Peter Clayson (3/1/15)
% peter.clayson@gmail.com

%see if the file be used for the figures and tables has been specified in
%varargin
if ~isempty(varargin)
    
    %check if data have been provided
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = varargin{ind+1};
        [pathpart,filepart] = fileparts(file);
    end
    
    try %see if ra_start is still open
        %if ra_start is open, close it
        if ishandle(varargin{3})
            close(varargin{3})
        end
    catch
    end

end %if ~isempty(varargin)

%if the file was not specified, prompt the user to indicate where the file
%is located.
if ~exist('file','var')
    [filepart, pathpart] = uigetfile({'*.mat','MAT-files (*.mat)'},'Data');

    if filepart == 0 
        errordlg('No file selected','File Error');
        era_start;
    end

    fprintf('\n\nLoading Data...\n\n');
end

%create a gui to allow the user to specify what aspects of the data will be
%viewed
era_startview_fig(filepart,pathpart);

end

function era_startview_fig(filepart,pathpart)
%Input
% filepart - filename to be loaded
% patherpart - path to the file
%
%Output
% No variables will be outputted to the Matlab workspace. Based on the
%  inputs from this gui, era_relfigures will be executed to display various
%  figures and tables (for more information about the tables and figures
%  see the user manual for the ERA toolbox)

%define parameters for figure position
figwidth = 550;
figheight = 550;

%define default inputs

%value to be used for dependability cutoff
inputs.depvalue = .70;

%for the inputs below a value of 1 (default) indicates that the
%figure/table should be viewed; a value of 0 indicates the figure/table
%should not be viewed

%figure that displays the dependability as the number of
%trials included in a subject's average increases
inputs.plotdep = 1;

%figure that shows the intraclass correlation coefficients
inputs.ploticc = 1;

%table displaying information about cutoffs based on dependability
inputs.inctrltable = 1;

%table displaying information about overall dependability with data
%including all trials
inputs.overalltable = 1;

%table displaying information about between- and within-person standard
%deviations
inputs.showstddevt = 1;

%figure depicting between-person standard deviations
inputs.showstddevf = 1;

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
uicontrol(era_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' filepart],...
    'Position',[0 row figwidth 25]);          

%next row
row = row - (rowspace*1.5);

%Print the text for dependability cutoff with a box for the user to specify
%the input
uicontrol(era_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Dependability Cutoff:',...
    'Position', [lcol row figwidth/4 25]);  

inputs.h(1) = uicontrol(era_gui,'Style','edit','fontsize',14,...
    'String',inputs.depvalue,...
    'Position', [rcol+5 row figwidth/4 25]);  

%next row
row = row - rowspace-10;

%indicate that a checked box means yes
uicontrol(era_gui,'Style','text','fontsize',12,...
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
uicontrol(era_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String','Would you like to plot Number of Trials v Dependability?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(2) = uicontrol(era_gui,'Style','checkbox',...
    'Value',inputs.plotdep,...
    'Position', [rcol row+12 figwidth/2 25]); 

%next row
row = row - rowspace;

%plot ICCs
uicontrol(era_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like to plot intraclass correlation coefficients?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(3) = uicontrol(era_gui,'Style','checkbox',...
    'Value',inputs.ploticc,...
    'Position', [rcol row+12 figwidth/2 25]); 

%next row
row = row - rowspace;

%dependability cutoff table
uicontrol(era_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like a table of event specific dependability information?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(4) = uicontrol(era_gui,'Style','checkbox',...
    'Value',inputs.inctrltable,...
    'Position', [rcol row+12 figwidth/2 25]); 

%next row
row = row - rowspace;

%overall dependability table
uicontrol(era_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like a table of overall dependability information?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(5) = uicontrol(era_gui,'Style','checkbox',...
    'Value',inputs.overalltable,...
    'Position', [rcol row+12 figwidth/2 25]); 

%next row
row = row - rowspace;

%between- and within-person standard deviation tables
str = ['Would you like a table of the overall relative '...
    'sizes of sources of variance?'];
uicontrol(era_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(6) = uicontrol(era_gui,'Style','checkbox',...
    'Value',inputs.showstddevt,...
    'Position', [rcol row+12 figwidth/2 25]);

%next row
row = row - rowspace;

%plot between-person standard deviations
str = ['Would you like a figure showing the between-person standard '...
    'deviations?'];
uicontrol(era_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    str,...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(7) = uicontrol(era_gui,'Style','checkbox',...
    'Value',inputs.showstddevf,...
    'Position', [rcol row+12 figwidth/2 25]);

%next row for buttons
row = row - rowspace*1.5;

%Create a back button that will take the user back to era_start
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 row figwidth/4 50],...
    'Callback',{@era_svb,era_gui}); 

%Create button that will check the inputs and begin processing the data
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [5*figwidth/8 row figwidth/4 50],...
    'Callback',{@era_svh,filepart,pathpart,inputs}); 

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
%Inputs
% varargin
%  3 - filename
%  4 - path to file
%  5 - inputs from gui

%need to take extension off file
filename = strsplit(varargin{3},'.');

%load the data to be viewed
REL = load(fullfile(varargin{4},[filename{1} '.mat']));
REL = REL.RELout;

%pull the inputs out of varargin
inputs = varargin{5};

%pass inputs from gui to era_relfigures
era_relfigures('data',REL,'depcutoff',str2double(inputs.h(1).String),...
    'plotdep',inputs.h(2).Value,'ploticc',inputs.h(3).Value,...
    'showinct',inputs.h(4).Value,'showoverallt',inputs.h(5).Value,...
    'showstddevt',inputs.h(6).Value,'plotbetstddev',inputs.h(7).Value);

end

