function overalltable = era_depoverallt(varargin)
%Display a table with the overall dependability information for data after
%applying the cutoff
%
%era_depoverallt('era_data',era_data,'gui',1);
%
%Last Modified 9/9/20
%
%Inputs
% era_data - ERA Toolbox data structure array. 
% gui - 0 for off, 1 for on
%
%Outputs
% overalltable - table displaying dependability information
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
% by Peter Clayson (7/24/16)
% peter.clayson@gmail.com
%
%1/19/17 PC
% updated copyright
%
%6/22/18 PC
% added standard deviation to table
%
%8/28/20 PC
% add subject-level reliability functionality
%
%9/9/20 PC
% changes for displaying difference score reliability

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
overalldep = {};
mintrl = {};
maxtrl = {};
meantrl = {};
medtrl = {};
stdtrl = {};
goodn = {};
badn = {};

if ~strcmp(era_data.rel.analysis,'ic_sserrvar')
    %put data together to display in tables
    for gloc=1:ngroups
        for eloc=1:nevents
            
            %label for group and/or event
            switch analysis
                case 1
                    label{end+1} = 'Measurement'; %#ok<*AGROW>
                case 2
                    label{end+1} = gnames{gloc};
                case 3
                    label{end+1} = enames{eloc};
                case 4
                    label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
            end
            
            %create a string with the dependability point estimate and credible
            %interval for overall data
            overalldep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                era_data.relsummary.group(gloc).event(eloc).dep.m,...
                era_data.relsummary.group(gloc).event(eloc).dep.ll,...
                era_data.relsummary.group(gloc).event(eloc).dep.ul);
            
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
        
        if strcmp(era_data.rel.analysis,'ic_diff')
        
        switch analysis
            case 3
                label{end+1} = 'diff score';
            case 4
                label{end+1} = [gnames{gloc} ' - diff score'];
        end
        
        mintrl{end+1} = '---';
        maxtrl{end+1} = '---';
        meantrl{end+1} = '---';
        medtrl{end+1} = '---';
        stdtrl{end+1} = '---';
        
        %pull good and bad ns
        goodn{end+1} = era_data.relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = '---';
        
        %create a string with the dependability point estimate and credible
        %interval for cutoff data
        overalldep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            era_data.relsummary.group(gloc).diffscore.pt,...
            era_data.relsummary.group(gloc).diffscore.ll,...
            era_data.relsummary.group(gloc).diffscore.ul);
    end
    end
    
    
    
elseif strcmp(era_data.rel.analysis,'ic_sserrvar')
    for gloc=1:ngroups
        for eloc=1:nevents
            
            bp_var = cell2mat(era_data.relsummary.data.g(gloc).e(eloc).gro_sds(:,1));
            pop_sdlog = era_data.relsummary.data.g(gloc).e(eloc).pop_sdlog;
            ciedge = .025;
            
            idtable_raw = era_data.relsummary.group(gloc).event(eloc).ssrel_table;
            
            goodid_rows = ismember(idtable_raw.id,...
                era_data.relsummary.group(gloc).goodids);
            
            idtable = idtable_raw(goodid_rows,:);
            
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
            
            %create a string with the dependability point estimate and credible
            %interval for overall data
            overalldep{end+1} = sprintf(' %0.2f SD: %0.2f',...
                mean(idtable.dep_pt),...
                std(idtable.dep_pt));
            
            %put together trial summary information
            mintrl{end+1} = min(idtable.trls);
            maxtrl{end+1} = max(idtable.trls);
            meantrl{end+1} = mean(idtable.trls);
            medtrl{end+1} = median(idtable.trls);
            stdtrl{end+1} = std(idtable.trls);
            
            %pull good and bad ns
            goodn{end+1} = sum(idtable.ind2include);
            
            badid_rows = ~ismember(idtable_raw.id,...
                era_data.relsummary.group(gloc).goodids);
            
            badn{end+1} = height(idtable_raw(badid_rows,:));
            
        end
    end
    
    
end
%create table to describe the data including all trials 
overalltable = table(label',goodn',badn',overalldep',meantrl',...
    medtrl',stdtrl',mintrl',maxtrl');

overalltable.Properties.VariableNames = {'Label', ...
    'n_Included','n_Excluded', ...
    'Dependability', 'Mean_Num_Trials', 'Med_Num_Trials',...
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
        'name','Dependability Analyses',...
        'numbertitle','off',...
        'resize','off');
    
    %Print the name of the loaded dataset
    uicontrol(era_overall,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String','Overall Dependability',...
        'Position',[0 row figwidth 25]);
    
    %Start a table
    t = uitable('Parent',era_overall,'Position',...
        [25 100 figwidth-50 figheight-175],...
        'Data',table2cell(overalltable));
    
    if ~strcmp(era_data.rel.analysis,'ic_sserrvar')
        set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
            'Dependability' 'Mean # of Trials' 'Med # of Trials'...
            'Std Dev of Trials' 'Min # of Trials' 'Max # of Trials'});
    elseif strcmp(era_data.rel.analysis,'ic_sserrvar')
        set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
            'Subject Level Dependability' 'Mean # of Trials' 'Med # of Trials'...
            'Std Dev of Trials' 'Min # of Trials' 'Max # of Trials'});
    end
    
    set(t,'ColumnWidth',{'auto' 'auto' 'auto' 110 'auto' 'auto' 'auto' 'auto'});
    set(t,'RowName',[]);
    
    %Create a save button that will take save the table
    uicontrol(era_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save Table',...
        'Position', [figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_saveoveralltable,era_data,overalltable});
    
    %Create button that will save good/bad ids
    uicontrol(era_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save IDs',...
        'Position', [5*figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_saveids,era_data});
    
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

%save either an excel or csv file
if strcmp(ext,'.xlsx')
    
    %print header information about the dataset
    filehead = {'Dependability Table Generated on'; datestr(clock);''}; 
    filehead{end+1} = sprintf('ERA Toolbox v%s',era_data.ver);
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',era_data.rel.filename);
    filehead{end+1} = sprintf('Dependability Cutoff: %0.2f',...
        era_data.relsummary.depcutoff);
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
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        era_data.relsummary.depcutoff);
    fprintf(fid,'Cutoff Threshold used the %s\n',...
        era_data.relsummary.meascutoff);
    fprintf(fid, 'Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s', strcat('Label,N Included,N Excluded,',...
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

%save the information in either an excel or csv format
if strcmp(ext,'.xlsx')
    
    datap{1,1} = 'Data Generated on';
    datap{end+1,1} = datestr(clock); 
    datap{end+1,1} = sprintf('ERA Toolbox v%s',era_data.ver);
    datap{end+1,1} = '';
    datap{end+1,1} = sprintf('Dataset: %s',era_data.rel.filename);
    datap{end+1,1} = sprintf('Dependability Cutoff: %0.2f',...
        era_data.relsummary.depcutoff);
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
    fprintf(fid,'Dependability Cutoff: %0.2f\n',...
        era_data.relsummary.depcutoff);
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
