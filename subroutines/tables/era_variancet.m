function vartable = era_variancet(varargin)
%Table of ICCs and between- and within-person standard deviations
%
%era_depvartable('era_data',era_data,'gui',1);
%
%Last Modified 8/3/19
%
%Inputs
% era_data - ERA Toolbox data structure array.
% gui - 0 for off, 1 for on
%
%Outputs
% vartable - table displaying information for ICCs and between- and
%  within-person standard deviations
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
%8/3/19 PC
% added button to view trt variance components

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
icc = {};
betsd = {};
witsd = {};

if ~strcmp(era_data.rel.analysis,'ic_sserrvar')
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
            
            %create a string with the icc point estimate and credible interval
            icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                era_data.relsummary.group(gloc).event(eloc).icc.m,...
                era_data.relsummary.group(gloc).event(eloc).icc.ll,...
                era_data.relsummary.group(gloc).event(eloc).icc.ul);
            
            %create a string with the between-person standard devation point
            %estimate and credible interval
            betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                era_data.relsummary.group(gloc).event(eloc).betsd.m,...
                era_data.relsummary.group(gloc).event(eloc).betsd.ll,...
                era_data.relsummary.group(gloc).event(eloc).betsd.ul);
            
            %create a string with the within-person standard devation point
            %estimate and credible interval
            witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                era_data.relsummary.group(gloc).event(eloc).witsd.m,...
                era_data.relsummary.group(gloc).event(eloc).witsd.ll,...
                era_data.relsummary.group(gloc).event(eloc).witsd.ul);
            
        end
    end
    
elseif strcmp(era_data.rel.analysis,'ic_sserrvar')
    
    %put data together to display in tables
    for gloc=1:ngroups
        for eloc=1:nevents
            
            bp_var = cell2mat(era_data.relsummary.data.g(gloc).e(eloc).gro_sds(:,1));
            pop_sdlog = era_data.relsummary.data.g(gloc).e(eloc).pop_sdlog;
            ciedge = .025;
            
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
            
            %create a string with the icc point estimate and credible interval
            icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                mean(bp_var ./ (bp_var + exp(pop_sdlog).^2)),...
                quantile(bp_var ./ (bp_var + exp(pop_sdlog).^2),ciedge),...
                quantile(bp_var ./ (bp_var + exp(pop_sdlog).^2),1-ciedge));
            
            %create a string with the between-person standard devation point
            %estimate and credible interval
            betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                sqrt(mean(bp_var)),...
                sqrt(quantile(bp_var,ciedge)),...
                sqrt(quantile(bp_var,1-ciedge)));
            
            %create a string with the within-person standard devation point
            %estimate and credible interval
            witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                exp(mean(pop_sdlog)),...
                exp(quantile(pop_sdlog,ciedge)),...
                exp(quantile(pop_sdlog,1-ciedge)));
            
        end
    end
    
    
end

%create table for displaying between-person and within-person standard
%deviation information
vartable = table(label',betsd',witsd',icc');

if ~strcmp(era_data.rel.analysis,'ic_sserrvar')
    vartable.Properties.VariableNames = {'Label',...
        'Between_StdDev','Within_StdDev','ICC'};
elseif strcmp(era_data.rel.analysis,'ic_sserrvar')
    vartable.Properties.VariableNames = {'Label',...
        'Group_Between_StdDev','Group_Within_StdDev',...
        'Group_ICC'};
end

%define parameters for figure size
figwidth = 674;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;
name = ['Point and 95% Interval Estimates for the Between-'...
    'and Within-Person Standard Deviations and ICCs'];

%create gui for standard-deviation table
var_gui= figure('unit','pix',...
    'position',[1250 600 figwidth figheight],...
    'menub','no',...
    'name',name,...
    'numbertitle','off',...
    'resize','off');

if ~strcmp(era_data.rel.analysis,'ic_sserrvar')
    %Print the name of the loaded dataset
    uicontrol(var_gui,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',...
        'Group-Level Between- and Within-Person Standard Deviations and ICCs',...
        'Position',[0 row figwidth 25]);
end

%Start a table
t = uitable('Parent',var_gui,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(vartable));
set(t,'ColumnName',{'Label' 'Between Std Dev'...
    'Within Std Dev' 'ICC'});
set(t,'ColumnWidth',{200 140 140 140});
set(t,'RowName',[]);

if ~strcmp(era_data.rel.analysis,'ic_sserrvar')
    
    %Create a save button that will save the table
    uicontrol(var_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save Table',...
        'Position', [figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_savevartable,era_data,vartable});
    
elseif strcmp(era_data.rel.analysis,'ic_sserrvar')
    
    %Create a save button that will save the table
    uicontrol(var_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','<html><tr><td align=center>Save Group-Level<br>Estimates',...
        'Position', [figwidth/8 25 figwidth/4 50],...
        'Callback',{@era_savevartable,era_data,vartable});
    
end

if strcmp(era_data.rel.analysis,'trt')
    
    %Create a save button that will view variance components
    uicontrol(var_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','<html><tr><td align=center>View All<br>Standard Deviations',...
        'Position', [(figwidth/8)*5 25 figwidth/4 50],...
        'Callback',{@era_viewvarcomp,era_data});
    
elseif strcmp(era_data.rel.analysis,'ic_sserrvar')
    
    %Create a save button that will view variance components
    uicontrol(var_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','<html><tr><td align=center>Save Subject<br>Coefficients',...
        'Position', [(figwidth/8)*5 25 figwidth/4 50],...
        'Callback',{@era_savesscoeffs,era_data});
    
end

end



function era_savevartable(varargin)
%if the user pressed the button to save the table with the information
%about between- and within-person standard deviations

%parse inputs
era_data = varargin{3};
vartable = varargin{4};

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
    filehead{end+1} = sprintf('Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(vartable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',era_data.ver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',era_data.rel.filename);
    fprintf(fid,'Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s', strcat('Label,Between-Person Std Dev',...
        ',Within-Person Std Dev,ICC'));
    fprintf(fid,' \n');
    
    for i = 1:height(vartable)
        formatspec = '%s,%s,%s,%s\n';
        fprintf(fid,formatspec,char(vartable{i,1}),...
            char(vartable{i,2}), char(vartable{i,3}),...
            char(vartable{i,4}));
    end
    
    fclose(fid);
    
end


end



function era_savesscoeffs(varargin)
%if the user pressed the button to save the table with the information
%about between- and within-person standard deviations

%parse inputs
era_data = varargin{3};

[savename, savepath] = uiputfile(...
    {'*.csv',...
    'Comma-Separated Vale File (.csv)'},...
    'Where would you like to save table?');

%save as csv file
fid = fopen(fullfile(savepath,savename),'w');
fprintf(fid,'%s\n','Table Generated on');
fprintf(fid,'%s\n',datestr(clock));
fprintf(fid,'ERA Toolbox v%s\n',era_data.ver);
fprintf(fid,' \n');
fprintf(fid,'Dataset: %s\n',era_data.rel.filename);
fprintf(fid,'Chains: %d, Iterations: %d\n',...
    era_data.rel.nchains,era_data.rel.niter);
fprintf(fid,'Reliability Cutoff: %0.2f\n',era_data.relsummary.depcutoff);
fprintf(fid,'%s\n','Subject-Specific Estimates');
fprintf(fid,' \n');
fprintf(fid,' \n');

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


if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

for gloc=1:ngroups
    for eloc=1:nevents
        %label for group and/or event
        switch analysis
            case 1
                label = 'Measurement';
            case 2
                label = gnames{gloc};
            case 3
                label = enames{eloc};
            case 4
                label = [gnames{gloc} ' - ' enames{eloc}];
        end
        
        t_out = era_data.relsummary.group(gloc).event(eloc).ssrel_table;
        
        fprintf(fid,'%s\n',label);
        
        fprintf(fid,'%s\n', strcat('ID,Participant Meets Reliability Cutoff',...
            ',Trial Counts,Dependability,Depenability Lower Limit',...
            ',Depenability Upper Limit,ICC,ICC Lower Limit,ICC Upper Limit',...
            ',Error Variance'));
        
        id = t_out.id;
        include = t_out.ind2include;
        trls = t_out.trls;
        dep_pt = t_out.dep_pt;
        dep_ll = t_out.dep_ll;
        dep_ul = t_out.dep_ul;
        icc_pt = t_out.icc_pt;
        icc_ll = t_out.icc_ll;
        icc_ul = t_out.icc_ul;
        ss_errvar = t_out.ss_errvar;
        
        
        for ii = 1:height(t_out)
            formatspec = '%s,%d,%d,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n';
            fprintf(fid,formatspec,...
                char(id{ii}),...
                include(ii),...
                trls(ii),...
                dep_pt(ii),...
                dep_ll(ii),...
                dep_ul(ii),...
                icc_pt(ii),...
                icc_ll(ii),...
                icc_ul(ii),...
                ss_errvar(ii));
        end
        
    end
end

fclose(fid);




end




function era_viewvarcomp(varargin)
%present another gui to view the variance components for trt data

%parse inputs
era_data = varargin{3};

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
bp = {};
bo = {};
bt = {};
txp = {};
oxp = {};
txo = {};
err = {};

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
        
        %create a string with the icc point estimate and credible interval
        bp{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_id.raw));
        
        bo{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw));
        
        bt{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw));
        
        txp{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw));
        
        oxp{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw));
        
        txo{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw));
        
        err{end+1} = sprintf('%0.2f',...
            mean(era_data.relsummary.data.g(gloc).e(eloc).sig_err.raw));
        
    end
end


%create table for displaying standard deviations
vartable = table(label',bp',bo',bt',oxp',txp',txo',err');

vartable.Properties.VariableNames = {'Label',...
    'bp','bo','bt','oxp','txp','txo','err'};

%define parameters for figure size
figwidth = 1000;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;
name = 'Point Estimates of All Standard Deviation Components';

%create gui for standard-deviation table
varall_gui= figure('unit','pix',...
    'position',[1250 600 figwidth figheight],...
    'menub','no',...
    'name',name,...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(varall_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    'Estimated Standard Deviation Components',...
    'Position',[0 row figwidth 25]);

%Start a table
t = uitable('Parent',varall_gui,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(vartable));
set(t,'ColumnName',{'Label'...
    'Between Person' 'Between Session' 'Between Trial' ...
    'Person x Session' 'Person x Trial' 'Session x Trial'...
    'Within Person (Error)'});
set(t,'ColumnWidth',{200 130 130 130 130 130 130 130});
set(t,'RowName',[]);

%Create a save button that will save the table
uicontrol(varall_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/8 25 figwidth/4 50],...
    'Callback',{@era_saveallvartable,era_data,vartable});

end

function era_saveallvartable(varargin)
%if the user pressed the button to save all of the standard deviation
%components

%parse inputs
era_data = varargin{3};
vartable = varargin{4};

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
    filehead{end+1} = sprintf('Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    filehead{end+1}='';
    filehead{end+1}='';
    
    xlswrite(fullfile(savepath,savename),filehead);
    writetable(vartable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'ERA Toolbox v%s\n',era_data.ver);
    fprintf(fid,' \n');
    fprintf(fid,'Dataset: %s\n',era_data.rel.filename);
    fprintf(fid,'Chains: %d, Iterations: %d',...
        era_data.rel.nchains,era_data.rel.niter);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    
    fprintf(fid,'%s', strcat('Label,Between Person, Between Session,',...
        'Between Trial,Person x Session,Person x Trial,Session x Trial,',...
        'Within Person (Error)'));
    fprintf(fid,' \n');
    
    for i = 1:height(vartable)
        formatspec = '%s,%s,%s,%s,%s,%s,%s,%s\n';
        fprintf(fid,formatspec,char(vartable{i,1}),...
            char(vartable{i,2}), char(vartable{i,3}),...
            char(vartable{i,4}), char(vartable{i,5}),...
            char(vartable{i,6}), char(vartable{i,7}),...
            char(vartable{i,8}));
    end
    
    fclose(fid);
    
end


end