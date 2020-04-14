function overalltable = era_trt_reloverallt(varargin)
%Display a table with the overall reliability information for data after
% applying the cutoff
%
%era_depoverallt('era_data',era_data,'gui',1);
%
%Last Modified 8/13/19
%
%Inputs
% era_data - ERA Toolbox data structure array.
% gui - 0 for off, 1 for on
%
%Outputs
% overalltable - table displaying reliability information
% a gui will also be shown if desired

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
% by Peter Clayson (8/2/19)
% peter.clayson@gmail.com
%
%8/13/19 PC
% Added option to specify the number of trials to use for computing retest
%  reliability
%
%8/20/19 PC
% Removed button to calculate retest reliability when viewing coefficients
%  of stability
%
%12/18/19 PC
% Fixed call to relcutoff instead of depcutoff in saving of files
%

%somersault through inputs
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help era_dep for more information about inputs'));
    end
    
    %check if era_data was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('era_data',varargin),1);
    if ~isempty(ind)
        era_data = varargin{ind+1};
    else
        error('varargin:era_data',... %Error code and associated error
            strcat('WARNING: era_data not specified \n\n',...
            'Please input era_data (ERA Toolbox data structure array).\n',...
            'See help era_depvtrialsplot for more information \n'));
    end
    
    %check if gui was specified
    %If it is not found, set display error.
    ind = find(strcmpi('gui',varargin),1);
    if ~isempty(ind)
        gui = varargin{ind+1};
    else
        error('varargin:gui',... %Error code and associated error
            strcat('WARNING: gui not specified \n\n',...
            'Please input gui specifying whether to display a gui.\n',...
            '0 for off, 1 for on\n',...
            'See help era_depvtrialsplot for more information \n'));
    end
    
end

%check whether any groups exist
if strcmpi(era_data.rel.groups,'none')
    ngroups = 1;
    %gnames = cellstr(era_data.rel.groups);
    gnames ={''};
else
    ngroups = length(era_data.rel.groups);
    gnames = era_data.rel.groups(:);
end

%check whether any events exist
if strcmpi(era_data.rel.events,'none')
    nevents = 1;
    %enames = cellstr(era_data.rel.events);
    enames = {''};
else
    nevents = length(era_data.rel.events);
    enames = era_data.rel.events(:);
end

%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

%create placeholders for displaying data in tables in guis
label = {};
overallrel = {};
mintrl = {};
maxtrl = {};
meantrl = {};
medtrl = {};
stdtrl = {};
goodn = {};
badn = {};

%put data together to display in tables
for gloc=1:ngroups
    for eloc=1:nevents
        
        %label for group and/or event
        switch analysis
            case 1
                label{end+1} = 'Measurement';
            case 2
                label{end+1} = gnames{gloc};
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end
        
        %create a string with the reliability point estimate and credible
        %interval for overall data
        overallrel{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            era_data.relsummary.group(gloc).event(eloc).rel.m,...
            era_data.relsummary.group(gloc).event(eloc).rel.ll,...
            era_data.relsummary.group(gloc).event(eloc).rel.ul);
        
        %put together trial summary information
        mintrl{end+1} = era_data.relsummary.group(gloc).event(eloc).trlinfo.min;
        maxtrl{end+1} = era_data.relsummary.group(gloc).event(eloc).trlinfo.max;
        meantrl{end+1} = era_data.relsummary.group(gloc).event(eloc).trlinfo.mean;
        medtrl{end+1} = era_data.relsummary.group(gloc).event(eloc).trlinfo.med;
        stdtrl{end+1} = era_data.relsummary.group(gloc).event(eloc).trlinfo.std;
        
        %pull good and bad ns
        goodn{end+1} = era_data.relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = length(era_data.relsummary.group(gloc).badids);
        
    end
end


%create table to describe the data including all trials
overalltable = table(label',goodn',badn',overallrel',meantrl',...
    medtrl',stdtrl',mintrl',maxtrl');

switch era_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

overalltable.Properties.VariableNames = {'Label', ...
    'n_Included','n_Excluded', ...
    rel_name, 'Mean_Num_Trials', 'Med_Num_Trials',...
    'Std_Num_Trials','Min_Num_Trials',...
    'Max_Num_Trials'};

%display gui if desired
if gui == 1
    
    %define parameters for figure size
    figwidth = 815;
    figheight = 500;
    
    %define space between rows and first row location
    rowspace = 25;
    row = figheight - rowspace*2;
    
    %create a gui for displaying the overall trial information
    era_overall= figure('unit','pix',...
        'position',[1150 150 figwidth figheight],...
        'menub','no',...
        'name',[rel_name ' Analyses Including All Trials'],...
        'numbertitle','off',...
        'resize','off');
    
    %Print the name of the loaded dataset
    uicontrol(era_overall,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',sprintf('Overall %s',rel_name),...
        'Position',[0 row figwidth 25]);
    
    %Start a table
    t = uitable('Parent',era_overall,'Position',...
        [25 100 figwidth-50 figheight-175],...
        'Data',table2cell(overalltable));
    set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
        rel_name 'Mean # of Trials' 'Med # of Trials'...
        'Std Dev of Trials' 'Min # of Trials' 'Max # of Trials'});
    set(t,'ColumnWidth',{'auto' 'auto' 'auto' 110 'auto' 'auto' 'auto' 'auto'});
    set(t,'RowName',[]);
    
    %Create a save button that will take save the table
    uicontrol(era_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save Table',...
        'Position', [.5*figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_saveoveralltable,era_data,overalltable});
    
    %Create button that will save good/bad ids
    uicontrol(era_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save IDs',...
        'Position', [3*figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_saveids,era_data});
    
    str = sprintf(['Use this button to estimate the test-retest reliability\n'...
        'for a given number of trials. This is helpful when adequate reliability\n'...
        'is not reached but you want to estimate obtained reliability']);
    
    
    if strcmp(era_data.relsummary.reltype_name,'trt')
        %Create button that will estimate a new trt reliability coefficient if
        %looking at retest data
        uicontrol(era_overall,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String','<html><tr><td align=center>Estimate New<br>Reliability Coefficient',...
            'Tooltip', str,...
            'Position', [5.5*figwidth/8 25 figwidth/4 50],...
            'Callback',{@era_newrel,era_data});
        
    end
end

end



function era_saveoveralltable(varargin)
%if the button to save the overall trial information table was pressed

%parse inputs
era_data = varargin{3};
overalltable = varargin{4};

%ask the user where the file should be saved
if ~ismac %macs can't use xlswrite
    [savename, savepath] = uiputfile(...
        {'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
else
    [savename, savepath] = uiputfile(...
        {'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save table?');
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

switch era_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

%save either an excel or csv file
if strcmp(ext,'.xlsx')
    
    %print header information about the dataset
    filehead = {'Dependability Table Generated on'; datestr(clock);''};
    filehead{end+1} = sprintf('ERA Toolbox v%s',era_data.ver);
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',era_data.rel.filename);
    filehead{end+1} = sprintf('%s Cutoff: %0.2f',...
        rel_name,...
        era_data.relsummary.relcutoff);
    filehead{end+1} = sprintf('Cutoff Threshold used the %s',...
        era_data.relsummary.meascutoff);
    filehead{end+1} = sprintf('Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    filehead{end+1}='';
    filehead{end+1}='';
    
    %write table
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(overalltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    %print header information about dataset
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s','Dependability Table Generated on ');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',era_data.ver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',era_data.rel.filename);
    fprintf(fid,'%s Cutoff: %0.2f\n',...
        rel_name,...
        era_data.relsummary.relcutoff);
    fprintf(fid,'Cutoff Threshold used the %s\n',...
        era_data.relsummary.meascutoff);
    fprintf(fid, 'Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s%s%s', strcat('Label,N Included,N Excluded,',...
        rel_name,...
        'Dependability,Mean Num of Trials,Med Num of Trials,',...
        'Std Dev Num of Trials,Min Num of Trials,Max Num of Trials'));
    fprintf(fid,' \n');
    
    %write the table information
    for i = 1:height(overalltable)
        formatspec = '%s,%d,%d,%s,%0.2f,%d,%0.2f,%d,%d\n';
        fprintf(fid,formatspec,char(overalltable{i,1}),...
            cell2mat(overalltable{i,2}),cell2mat(overalltable{i,3}),...
            char(overalltable{i,4}),cell2mat(overalltable{i,5}),...
            cell2mat(overalltable{i,6}),cell2mat(overalltable{i,7}),...
            cell2mat(overalltable{i,8}),cell2mat(overalltable{i,9}));
    end
    
    fclose(fid);
    
end

end

function era_saveids(varargin)
%if the user pressed the button to save which ids were considered good and
%which were considered bad (based on whether data met the cutoff thresholds

era_data = varargin{3};

%ask the user where the file should be saved
if ~ismac
    [savename, savepath] = uiputfile(...
        {'*.xlsx','Excel File (.xlsx)';'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save the data?');
else
    [savename, savepath] = uiputfile(...
        {'*.csv',...
        'Comma-Separated Vale File (.csv)'},...
        'Where would you like to save the data?');
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

switch era_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

%save the information in either an excel or csv format
if strcmp(ext,'.xlsx')
    
    datap{1,1} = 'Data Generated on';
    datap{end+1,1} = datestr(clock);
    datap{end+1,1} = sprintf('ERA Toolbox v%s',era_data.ver);
    datap{end+1,1} = '';
    datap{end+1,1} = sprintf('Dataset: %s',era_data.rel.filename);
    datap{end+1,1} = sprintf('%s Cutoff: %0.2f',...
        rel_name,...
        era_data.relsummary.relcutoff);
    datap{end+1,1} = sprintf('Cutoff Threshold used the %s',...
        era_data.relsummary.meascutoff);
    datap{end+1,1} = sprintf('Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    datap{end+1,1}='';
    datap{end+1,1}='';
    datap{end+1,1} = 'Good IDs'; datap{end,2} = 'Bad IDs';
    
    gids = [];
    bids = [];
    srow = length(datap);
    
    for j=1:length(era_data.relsummary.group)
        gids = [gids;era_data.relsummary.group(j).goodids(:)];
        bids = [bids;era_data.relsummary.group(j).badids(:)];
    end
    
    for i = 1:length(gids)
        datap{i+srow,1}=char(gids(i));
    end
    
    for i = 1:length(bids)
        datap{i+srow,2}=char(bids(i));
    end
    
    xlswrite(fullfile(savepath,savename),datap);
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Data Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',era_data.ver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',era_data.rel.filename);
    fprintf(fid,'%s Cutoff: %0.2f\n',...
        rel_name,...
        era_data.relsummary.relcutoff);
    fprintf(fid,'Cutoff Threshold used the %s\n',...
        era_data.relsummary.meascutoff);
    fprintf(fid,'Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    fprintf(fid,'%s\n','Good IDs,Bad IDs');
    
    gids = [];
    bids = [];
    
    for j=1:length(era_data.relsummary.group)
        gids = [gids;era_data.relsummary.group(j).goodids(:)];
        bids = [bids;era_data.relsummary.group(j).badids(:)];
    end
    
    maxlength = max([length(gids) length(bids)]);
    minlength = min([length(gids) length(bids)]);
    
    if minlength == length(gids)
        whichlonger = 1;
    elseif maxlength == length(gids)
        whichlonger = 2;
    end
    
    for i = 1:maxlength
        if i <= minlength
            fprintf(fid,'%s,%s\n',gids{i},bids{i});
        elseif i > minlength && whichlonger == 1
            fprintf(fid,',%s\n',bids{i});
        elseif i > minlength && whichlonger == 2
            fprintf(fid,'%s,\n',gids{i});
        end
    end
    
    fclose(fid);
    
end


end



function era_newrel(varargin)
%if the user pressed the button to calculate a new test-retest reliability
%coefficient, this gui will be pulled up

era_data = varargin{3};

%define parameters for figure position
figwidth = 550;
figheight = 250;

%define space between rows and first row location
rowspace = 40;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

%specify font size
fsize = 14;

%create the gui
era_newrelgui = figure('unit','pix',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'name','Specify Inputs',...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(era_newrelgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' era_data.rel.filename],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

%next row
row = row - (rowspace*.45);

%Print the name of the measurement analyzed
uicontrol(era_newrelgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',['Measurement:  ' era_data.proc.measheader],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

%next row
row = row - (rowspace*1.4);

str = sprintf(['Number of trials to use for recalculating reliability\n',...
    'This should be the mean or median number of trials retained for averaging']);

%Print the text for reliability cutoff with a box for the user to specify
%the input
uicontrol(era_newrelgui,...
    'Style','text',...
    'fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Number of Trials:',...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);

inputs.trls = uicontrol(era_newrelgui,...
    'Style','edit',...
    'fontsize',fsize,...
    'String','',...
    'Position', [rcol+5 row+6 figwidth/4 25]);

%next row
row = row - rowspace*2;

str = sprintf(['Use this button to estimate the new test-retest reliability\n'...
    'information for the specified number of trials.']);

%Create button that will save good/bad ids
uicontrol(era_newrelgui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','<html><tr><td align=center>Estimate New<br>Relability Coefficient',...
    'Tooltip', str,...
    'Position', [.5*figwidth/8 row figwidth/1.11 50],...
    'Callback',{@era_shownewrel,era_data,inputs.trls});

end

function era_shownewrel(varargin)
%estimate the new reliability estimates and show it

era_data = varargin{3};
ntrls = str2double(varargin{4}.String);

%place era_data.data in REL to work with
REL = era_data.rel;

%check whether any groups exist
if strcmpi(REL.groups,'none')
    ngroups = 1;
    %gnames = cellstr(REL.groups);
    gnames ={''};
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

%check whether any events exist
if strcmpi(REL.events,'none')
    nevents = 1;
    %enames = cellstr(REL.events);
    enames = {''};
else
    nevents = length(REL.events);
    enames = REL.events(:);
end

%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider


if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

%extract information from REL and store in data for crunching
switch analysis
    case 1 %1 - no groups or event types to consider
        
        gloc = 1;
        eloc = 1;
        
        data.g(gloc).e(eloc).label = REL.out.labels{gloc};
        data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
        data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
        data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
        data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
        data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
        data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
        data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
        data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
        data.g(gloc).e(eloc).elabel = cellstr('none');
        data.g(gloc).glabel = gnames(gloc);
        
    case 2 %2 - possible multiple groups but no event types to consider
        
        eloc = 1;
        
        for gloc=1:length(REL.out.labels)
            
            data.g(gloc).e(eloc).label = REL.out.labels{gloc};
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
            data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
            data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
            data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
            data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
            data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
            data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
            data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
            data.g(gloc).e(eloc).elabel = cellstr('none');
            data.g(gloc).glabel = gnames(gloc);
            
        end
        
    case 3 %3 - possible event types but no groups to consider
        
        gloc = 1;
        
        for eloc=1:length(REL.out.labels)
            
            data.g(gloc).e(eloc).label = REL.out.labels{eloc};
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
            data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,eloc);
            data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,eloc);
            data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,eloc);
            data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,eloc);
            data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,eloc);
            data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,eloc);
            data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,eloc);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
            
        end
        
    case 4 %4 - possible groups and event types to consider
        for ii=1:length(REL.out.labels)
            
            %use the underscores that were added in era_computerel to
            %differentiate where the group and event the data are for
            lblstr = strsplit(REL.out.labels{ii},'_;_');
            
            eloc = find(ismember(enames,lblstr(2)));
            gloc = find(ismember(gnames,lblstr(1)));
            
            data.g(gloc).e(eloc).label = REL.out.labels{ii};
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,ii);
            data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,ii);
            data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,ii);
            data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,ii);
            data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,ii);
            data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,ii);
            data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,ii);
            data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,ii);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
            
        end
end %switch analysis

%compute reliability data for each group and event
switch analysis
    case 1 %no groups or event types to consider
        
        %the same generic structure is used for relsummary, so the event
        %and group locations will both be 1 for the data
        eloc = 1;
        gloc = 1;
        
        %compute reliabiltiy
        [llrel,mrel,ulrel] = era_rel_trt(...
            'gcoeff',era_data.relsummary.gcoeff,...
            'reltype',era_data.relsummary.reltype,...
            'bp',data.g(gloc).e(eloc).sig_id.raw,...
            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
            'err',data.g(gloc).e(eloc).sig_err.raw,...
            'obs', ntrls,'CI',era_data.relsummary.ciperc);
        
        newrelsummary.group(gloc).event(eloc).rel.m = mrel;
        newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
        newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
        
    case 2 %possible multiple groups but no event types to consider
        
        
        %since the same generic structure is used for relsummary, the event
        %location will be defined as 1.
        eloc = 1;
        
        for gloc=1:ngroups %loop through each group
            
            %compute reliabiltiy
            [llrel,mrel,ulrel] = era_rel_trt(...
                'gcoeff',era_data.relsummary.gcoeff,...
                'reltype',era_data.relsummary.reltype,...
                'bp',data.g(gloc).e(eloc).sig_id.raw,...
                'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                'err',data.g(gloc).e(eloc).sig_err.raw,...
                'obs',ntrls,'CI',era_data.relsummary.ciperc);
            
            newrelsummary.group(gloc).event(eloc).rel.m = mrel;
            newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
            newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
            
        end
        
    case 3 %possible event types but no groups to consider
        
        %since the relsummary structure is generic for any number of groups
        %or events, the data will count as being from 1 group.
        gloc = 1;
        
        %cylce through each event
        for eloc=1:nevents
            
            %compute reliabiltiy
            [llrel,mrel,ulrel] = era_rel_trt(...
                'gcoeff',era_data.relsummary.gcoeff,...
                'reltype',era_data.relsummary.reltype,...
                'bp',data.g(gloc).e(eloc).sig_id.raw,...
                'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                'err',data.g(gloc).e(eloc).sig_err.raw,...
                'obs',ntrls,'CI',era_data.relsummary.ciperc);
            
            newrelsummary.group(gloc).event(eloc).rel.m = mrel;
            newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
            newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
            
        end
        
    case 4 %groups and event types to consider
        
        %cycle through each event
        for eloc=1:nevents
            
            %cycle through each group
            for gloc=1:ngroups
                
                %compute reliability
                [llrel,mrel,ulrel] = era_rel_trt(...
                    'gcoeff',era_data.relsummary.gcoeff,...
                    'reltype',era_data.relsummary.reltype,...
                    'bp',data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',ntrls,'CI',era_data.relsummary.ciperc);
                
                newrelsummary.group(gloc).event(eloc).rel.m = mrel;
                newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
                newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
                
            end
        end
        
        
end %switch analysis


%create placeholders for displaying data in tables in guis
label = {};
overallrel = {};
new_trls = {};

%put data together to display in tables
for gloc=1:ngroups
    for eloc=1:nevents
        
        %label for group and/or event
        switch analysis
            case 1
                label{end+1} = 'Measurement';
            case 2
                label{end+1} = gnames{gloc};
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end
        
        %create a string with the reliability point estimate and credible
        %interval for overall data
        overallrel{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            newrelsummary.group(gloc).event(eloc).rel.m,...
            newrelsummary.group(gloc).event(eloc).rel.ll,...
            newrelsummary.group(gloc).event(eloc).rel.ul);
        
        new_trls{end+1} = ntrls;
        
    end
end


%create table to describe the data including all trials
overalltable = table(label',overallrel',new_trls');

switch era_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

overalltable.Properties.VariableNames = {'Label', ...
    rel_name, 'Num_Trials'};

%define parameters for figure size
figwidth = 415;
figheight = 350;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%create a gui for displaying the overall trial information
era_overall= figure('unit','pix',...
    'position',[1150 150 figwidth figheight],...
    'menub','no',...
    'name',sprintf('%s Analyses for %d Trials',...
    rel_name,ntrls),...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(era_overall,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',sprintf('Overall %s',rel_name),...
    'Position',[0 row figwidth 25]);

%Start a table
t = uitable('Parent',era_overall,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(overalltable));
set(t,'ColumnName',{'Label' ...
    rel_name '# of Trials'});
set(t,'ColumnWidth',{151 110 100});
set(t,'RowName',[]);


end