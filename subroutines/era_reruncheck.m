function era_reruncheck
%
%Execute gui to ask user whether to rerun model with more iterations
%
%Last Updated 8/16/17
%
%
%Input
% No inputs required
%
%Output
% No output generated. The user's choice will be checked in era_startproc
%

% Copyright (C) 2016-2017 Peter E. Clayson
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
%1/19/17 PC
% updated copyright
%
%8/16/17 PC
% fixed typo in comments

%define parameters for figure position
figwidth = 500;
figheight = 200;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2.5;

%initialize gui
era_gui= figure('unit','pix','Visible','off',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'numbertitle','off',...
  'resize','off');

movegui(era_gui,'center');

str = {'Chains did not converge';...
    'Would you like to rerun with more iterations?'};

%Write text
uicontrol(era_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',str,...
    'Position',[0 row figwidth 50]);          

%Create a button that will not rerun the model
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','<html><center>Do Not<br>Rerun',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_setrerun,0,era_gui}); 

%Create button that will rerun the model with more iterations
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Rerun',...
    'Position', [5*figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_setrerun,1,era_gui}); 

%display gui
set(era_gui,'Visible','on');

%tag the gui
era_gui.Tag = 'era_gui';

%wait so the user can respond
uiwait(era_gui);


end

function era_setrerun(varargin)
%to rerun or not to rerun the model... that is the question :)
uiresume;
era_gui = findobj('Tag','era_gui');
guidata(era_gui,varargin{3});

end


