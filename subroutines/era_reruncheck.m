function era_reruncheck
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

%Create a button that will take the user to the gui for setting the inputs
%to process data
uicontrol(era_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','<html><center>Do Not<br>Rerun',...
    'Position', [figwidth/8 25 figwidth/3 75],...
    'Callback',{@era_setrerun,0,era_gui}); 

%Create button that will take the user to the gui for setting the inputs
%for viewing the data
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

uiresume;
era_gui = findobj('Tag','era_gui');
guidata(era_gui,varargin{3});

end


