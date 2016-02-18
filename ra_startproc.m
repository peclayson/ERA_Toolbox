function ra_startproc

[filepart, pathpart] = uigetfile('*.*',...
    'Select the data file to be processed');

if filepart == 0 
    errordlg('No file selected','File Error');
    ra_start;
end

fprintf('\n\nLoading Data...\n\n');

%load file
dataraw = readtable(fullfile(pathpart,filepart));

%pull the headernames from the file
collist = dataraw.Properties.VariableNames;
collist{end+1} = 'none';

ra_startproc_fig(collist,filepart,pathpart,dataraw);

end

function ra_startproc_fig(collist,filepart,pathpart,dataraw)

%define parameters for figure position
figwidth = 500;
figheight = 300;
collist_nonone = collist;
collist_nonone(end) = [];

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/2);

ra_gui= figure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Specify Inputs for Analysis',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(ra_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' filepart],...
    'Position',[0 row figwidth 25]);          

%next row
row = row - rowspace*1.5;

%Print the headers
uicontrol(ra_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Variable',...
    'Position', [figwidth/8 row figwidth/4 25]);  

uicontrol(ra_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Input',...
    'Position',[5*(figwidth/8) row figwidth/4 25]);

%next row
row = row - rowspace*1.5;

%Print the name of the variables and place a listbox with possible options
uicontrol(ra_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Participant ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(1) = uicontrol(ra_gui,'Style','pop','fontsize',14,...
    'String',collist_nonone,'Value',1,...
    'Position', [rcol row figwidth/2 25]);  

row = row - rowspace;

uicontrol(ra_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Measurement:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(2) = uicontrol(ra_gui,'Style','pop','fontsize',14,...
    'String',collist_nonone,'Value',2,...
    'Position', [rcol row figwidth/2 25]); 

row = row - rowspace;

uicontrol(ra_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Group ID:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(3) = uicontrol(ra_gui,'Style','pop','fontsize',14,...
    'String',collist,'Value',length(collist),...
    'Position', [rcol row figwidth/2 25]); 

row = row - rowspace;

uicontrol(ra_gui,'Style','text','fontsize',14,...
    'HorizontalAlignment','left',...
    'String','Event Type:',...
    'Position', [lcol row figwidth/4 25]);  

inplists(4) = uicontrol(ra_gui,'Style','pop','fontsize',14,...
    'String',collist,'Value',length(collist),...
    'Position', [rcol row figwidth/2 25]); 


row = row - rowspace*2.5;

%Create a back button that will take the user back to ra_start
uicontrol(ra_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Back',...
    'Position', [figwidth/8 row figwidth/4 50],...
    'Callback',@bb_call); 

%Create button that will check the inputs and begin processing the data
uicontrol(ra_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Analyze',...
    'Position', [5*figwidth/8 row figwidth/4 50],...
    'Callback',{@ra_exec,inplists,collist,filepart,pathpart,dataraw}); 


end

%if back button is pushed, go back to ra_start
function bb_call 
   
ra_start;
   
end

%if execute button is pushed, analyze the loaded data
function ra_exec(varargin)

inp = varargin{3};
collist = varargin{4};
filepart = varargin{5};
pathpart = varargin{6};
dataraw = varargin{7};

choices = cell2mat(get(inp(:),'value'));

%check if there are duplicates (other than 'none')
[n, bin] = histc(choices, unique(choices));
multiple = find(n > 1);
ind = find(ismember(bin, multiple));
if ~isempty(ind)
    probcol = {};
    for i = 1:length(ind)
        if choices(ind(i)) ~= length(collist)
            probcol(end+1) = collist(ind(i));
        end
    end
    if ~isempty(probcol)
        dlg = {'Duplicate variable names were not provided for ';...
            probcol; ...
            'When selecting column headers, please select unique names'};
        errordlg(dlg, 'Unique variable names not provided');
        ra_startproc_fig(collist,filepart);
    end
end

%parse inputs
idheader = char(collist(choices(1)));
measheader = char(collist(choices(2)));

if choices(3) ~= length(collist)
    groupheader = char(collist(choices(3)));
elseif choices(3) == length(collist)
    groupheader = '';
end

if choices(4) ~= length(collist)
    eventheader = char(collist(choices(4)));
elseif choices(4) == length(collist)
    eventheader = '';
end

closereq;

[savename, savepath] = uiputfile(fullfile(pathpart,'*.mat'),...
    'Where would you like to save the output files?');

dataout = ra_loadfile('file',fullfile(pathpart,filepart),...
    'idcol',idheader,'meascol',measheader,'groupcol',groupheader,...
    'eventcol',eventheader,'dataraw',dataraw);

%Change working dir for Stan files

mkdir(savepath,'Temp_StanFiles');

origdir = cd(fullfile(savepath,'Temp_StanFiles'));

RELout = ra_computerel('data',dataout);

cd(origdir);
rmdir(fullfile(savepath,'Temp_StanFiles'),'s');

fprintf('\n\nSaving Processed Data...\n\n');

%Matlab spits out various warnings when trying to save the StanFit part of
%the structure. Turning off warnings brielfy so the user does not become
%concerned
warning('off','all');
save(fullfile(savepath,savename),'RELout');
warning('on','all');

ra_startview(RELout);

end


