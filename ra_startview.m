function ra_startview(varargin)

if ~isempty(varargin)
    
    %check if data have been provided
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = varargin{ind+1};
        [pathpart,filepart] = fileparts(file);
    end
 
end %if ~isempty(varargin)

if ~exist('file','var')
    [filepart, pathpart] = uigetfile({'*.mat','MAT-files (*.mat)'},'Data');

    if filepart == 0 
        errordlg('No file selected','File Error');
        ra_start;
    end

    fprintf('\n\nLoading Data...\n\n');
end

ra_startview_fig(filepart,pathpart);

end

function ra_startview_fig(filepart,pathpart)

%define parameters for figure position
figwidth = 550;
figheight = 450;

%default inputs
inputs.depvalue = .70;
inputs.plotdep = 1;
inputs.ploticc = 1;
inputs.inctrltable = 1;
inputs.overalltable = 1;

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

ra_gui= figure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Inputs',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(ra_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' filepart],...
    'Position',[0 row figwidth 25]);          

%next row
row = row - (rowspace*1.5);

%Print the name of the variables and place a listbox with possible options
uicontrol(ra_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Dependability Cutoff:',...
    'Position', [lcol row figwidth/4 25]);  

inputs.h(1) = uicontrol(ra_gui,'Style','edit','fontsize',14,...
    'String',inputs.depvalue,...
    'Position', [rcol+5 row figwidth/4 25]);  

row = row - rowspace-10;

uicontrol(ra_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','center',...
    'String','Checked = YES',...
    'Position', [rcol+5 row figwidth/4 25]);  

row = row - rowspace;
rowspace = 50;
rcol = (figwidth/4)*3;

uicontrol(ra_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String','Would you like to plot Number of Trials v Dependability?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(2) = uicontrol(ra_gui,'Style','checkbox',...
    'Value',inputs.plotdep,...
    'Position', [rcol row+12 figwidth/2 25]); 

row = row - rowspace;

uicontrol(ra_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like to plot intraclass correlation coefficients?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(3) = uicontrol(ra_gui,'Style','checkbox',...
    'Value',inputs.ploticc,...
    'Position', [rcol row+12 figwidth/2 25]); 

row = row - rowspace;

uicontrol(ra_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like a table of event specific dependability information?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(4) = uicontrol(ra_gui,'Style','checkbox',...
    'Value',inputs.inctrltable,...
    'Position', [rcol row+12 figwidth/2 25]); 

row = row - rowspace;

uicontrol(ra_gui,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',...
    'Would you like a table of overall dependability information?',...
    'Position', [lcol row figwidth/2 40]);  

inputs.h(5) = uicontrol(ra_gui,'Style','checkbox',...
    'Value',inputs.overalltable,...
    'Position', [rcol row+12 figwidth/2 25]); 

row = row - rowspace*1.5;

%Create a back button that will take the user back to ra_start
uicontrol(ra_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 row figwidth/4 50],...
    'Callback',{@ra_svb,ra_gui}); 

%Create button that will check the inputs and begin processing the data
uicontrol(ra_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [5*figwidth/8 row figwidth/4 50],...
    'Callback',{@ra_svh,filepart,pathpart,inputs,ra_gui}); 

end

function ra_svb(varargin)

close(varargin{3});
ra_start;

end

function ra_svh(varargin)

filename = strsplit(varargin{3},'.');

REL = load(fullfile(varargin{4},[filename{1} '.mat']));
REL = REL.RELout;

inputs = varargin{5};


ra_relfigures('data',REL,'relcutoff',str2double(inputs.h(1).String),...
    'plotdep',inputs.h(2).Value,'ploticc',inputs.h(3).Value,...
    'showinct',inputs.h(4).Value,'showoverallt',inputs.h(5).Value);

end

