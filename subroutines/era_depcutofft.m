function inctrltable = era_depcutofft(varargin)
%Display a table with the information for the cutoffs stratified by group
% and condition
%
%era_depcutofft('era_data',era_data,'gui',1);
%
%Last Modified 1/19/17
%
%Inputs
% era_data - ERA Toolbox data structure array. 
% gui - 0 for off, 1 for on
%
%Outputs
% inctrltable - table displaying the trial cutoff information
% a gui will also be shown if desired

% Copyright (C) 2016-2019 Peter E. Clayson
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
trlcutoff = {};
dep = {};

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
        
        %pull the trial cutoff
        trlcutoff{end+1} = era_data.relsummary.group(gloc).event(eloc).trlcutoff;
        
        %create a string with the dependability point estimate and credible
        %interval for cutoff data
        dep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            era_data.relsummary.group(gloc).event(eloc).rel.m,...
            era_data.relsummary.group(gloc).event(eloc).rel.ll,...
            era_data.relsummary.group(gloc).event(eloc).rel.ul);
        
    end 
end

%create table to display the cutoff information with its associated
%dependability info
inctrltable = table(label',trlcutoff',dep');

inctrltable.Properties.VariableNames = {'Label','Trial_Cutoff',...
    'Dependability'};

%only display gui if desired
if gui == 1

    %define parameters for figure size
    figwidth = 500;
    figheight = 400;

    %define space between rows and first row location
    rowspace = 25;
    row = figheight - rowspace*1.6;

    %create a gui to display the cutoff information table
    era_inctrl= figure('unit','pix',...
      'position',[1150 700 figwidth figheight],...
      'menub','no',...
      'name',...
      'Results of Increasing Trials the Number of Trials on Dependability',...
      'numbertitle','off',...
      'resize','off');

    %Add table title
    uicontrol(era_inctrl,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',...
        sprintf('Dependability Analyses, %0.2f Cutoff\n',...
        era_data.relsummary.depcutoff),...
        'Position',[0 row figwidth 20]);   

    uicontrol(era_inctrl,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',...
        sprintf('Cutoff Used the %s',...
        era_data.relsummary.meascutoff),...
        'Position',[0 row-20 figwidth 20]); 

    %Start a table
    t = uitable('Parent',era_inctrl,'Position',...
        [25 100 figwidth-50 figheight-175],...
        'Data',table2cell(inctrltable));
    set(t,'ColumnName',{'Label' 'Trial Cutoff' 'Dependability'});
    set(t,'ColumnWidth',{200 'auto' 170});
    set(t,'RowName',[]);

    %Create a save button that will take save the table
    uicontrol(era_inctrl,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save Table',...
        'Position', [figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_saveinctrltable,era_data,inctrltable}); 
end

end


function era_saveinctrltable(varargin)
%if the user pressed the button to save the table with the information
%about the cutoffs

%parse inputs
era_data = varargin{3};
incltrltable = varargin{4};

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

%save as either excel or csv file
if strcmp(ext,'.xlsx')
    
    filehead = {'Table Generated on'; datestr(clock);''}; 
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
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(incltrltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
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
    
    fprintf(fid,'%s', strcat('Label,Trial Cutoff,Dependability'));
    fprintf(fid,' \n');
    
    for i = 1:height(incltrltable)
         formatspec = '%s,%d,%s\n';
         fprintf(fid,formatspec,char(incltrltable{i,1}),...
             cell2mat(incltrltable{i,2}),char(incltrltable{i,3}));
    end
    
    fclose(fid);
    
end


end

