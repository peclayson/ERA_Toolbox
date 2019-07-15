function RELout = era_computerel(varargin)
%Prepare and execute cmdstan code for dependability analyses
%
%era_computerel('data',era_datatable,'chains',3,'iter',1000)
%
%Lasted Updated 4/18/19
%
%Required Inputs:
% data - data table outputted from the era_loadfile script (see era_loadfile
%  for more information about table format)
%  Note: era_loadfile sets up the datatable in a specific format with
%  specific header names that are used in this script
% chains - number of chains to run in stan 
% iter - number of iterations to run in stan 
%
%Optional Inputs:
% verbose - 1: Do not print iterations, 2: Print iterations (default: 1)
% showgui - 1: Do not show a gui while computations are running, 2: show a 
%  gui while the computations are running (default: 1)
%
%Outputs:
% RELout - structure array with the following fields.
%  filename: filename of the processed dataset
%  niter: number of interations run in Stan
%  events: names of the events (if applicable)
%  groups: names of the groups (if applicable)
%  data: compiled tables that were used to create data structures to pass
%   to stan
%  out: data output from stan and labels corresponding to order of data
%   mu: person mean
%   sig_u: person variance
%   sig_e: error variance
%   labels: labels for the data. If either events or groups were processed
%    then this the labeles for events or groups. If both events and groups
%    were processed then this labels will be 'event_group' for each
%    possible combination (e.g., 'Error_;_Controls' or 'Hits_;_Patients'  
%    event and group will be separated by _;_ to avoid using a common 
%    delimeter that could be used by a researcher to label an event)
%  conv: nest for convergence data
%   data: convergence statistics with three columns (parameter, n_eff,
%    r_hat)
%  stan_in: compiled model syntax that was passed to Stan

%For more information about how variance components were estimated please
% see
%
% Clayson, P. E., & Miller, G. A. (2017). ERP Reliability Analysis
% (ERA) Toolbox: An open-source toolbox for analyzing the reliability of
% event-related potentials. International Journal of Psychophysiology, 111,
% 68-79. doi: 10.1016/j.ijpsycho.2016.10.012
%
% Baldwin, S. A., Larson, M. J., & Clayson, P. E. (2015). The dependability 
% of electrophysiological measurements of performance monitoring in a 
% clinical sample: A generalizability and decision analysis of the ERN and 
% Pe. Psychophysiology, 52(6), 790-800. http://doi.org/10.1111/psyp.12401
%
%It would be terribly remiss of me not thank Dr. Scott Baldwin for
%conceptualizing and developing the formulas that are implemented in this
%toolbox. Dr. Baldwin also wrote the original Stan syntax in R and 
%graciously provided me with all of his code. This Matlab code is based
%off of his R code. 

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
% by Peter Clayson (4/18/16)
% peter.clayson@gmail.com
%
%7/19/16 PC
% Changes associate with updating to cmdstan-2.10.0
%  changed <- to = (<- is not deprecated)
%
%7/28/16 PC
% Added defineversion function
%
%7/30/16 PC
% Removed some old code that is no longer used
% Add the option to pop-up a gui while code is running
%
%10/22/16 PC
% Vectorized code to improve speed
%
%1/19/17 PC
% updated copyright
%
%6/25/17 PC
% added error check for missing cells in dataset
%
%2/28/19 PC
% added functionality for computing test-retest reliability
%
%4/18/19 PC
% computing all facets simultaneously for trt analyses took much too long,
%  so I broke it up by task and event which cut down processing
%  substantially (from ~36 hours to ~12 hours)
%
%6/17/19 PC
% added seed to all cmdstan analyses for reproducibility of results

%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_computerel for more information on optional inputs'));
    end
    
    %check if the dataset is present
    ind = find(strcmp('data',varargin),1);
    if ~isempty(ind)
        datatable = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
    
    %check whether chains is specified
    ind = find(strcmp('chains',varargin),1);
    if ~isempty(ind)
        nchains = varargin{ind+1}; 
    else 
        error('varargin:nchains',... %Error code and associated error
        strcat('WARNING: Number of chains not specified \n\n',... 
        'Please input the number of chains for stan\n',...
        'See help era_computerel for more information on inputs'));
    end
    
    %check whether iter is specified
    ind = find(strcmp('iter',varargin),1);
    if ~isempty(ind)
        niter = varargin{ind+1}; 
    else 
        error('varargin:niter',... %Error code and associated error
        strcat('WARNING: Number of iterations not specified \n\n',... 
        'Please input the number of iterations for stan\n',...
        'See help era_computerel for more information on inputs'));
    end
    
    %check whether verbose is specified
    ind = find(strcmp('verbose',varargin),1);
    if ~isempty(ind)
        verbose = varargin{ind+1}; 
    else 
        verbose = 1;
    end
    
    %check whether showgui is specified
    ind = find(strcmp('showgui',varargin),1);
    if ~isempty(ind)
        showgui = varargin{ind+1}; 
    else 
        showgui = 1;
    end
    
elseif isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help era_computerel for more information on inputs'));
    
end %if ~isempty(varargin)

eraver = era_defineversion;

%ensure the necessary columns are present in the table (at least the
%headers for id and meas)
colnames = datatable.Properties.VariableNames;

%check for id column
if ~sum(strcmpi(colnames,'id')) 
    if ~exist('headererror','var')
        headerror{1} = 'Subject ID';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Subject ID';
    end
end

%check for measurement column
if ~sum(strcmpi(colnames,'meas')) 
    if ~exist('headererror','var')
        headerror{1} = 'Measurement';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Measurement';
    end
end

%error catch in case headers for the columns needed are not specified
if exist('headerror','var')
    error('varargin:colheaders',... %Error code and associated error
    strcat('WARNING: Column headers not properly defined in data \n\n',... 
    'Please properly specify the headers for\n',...
    char(strjoin(headerror,', ')),'\n',...
    'See help era_loadfile for data table format \n'));
end

%check whether there are missing cells in the datatable. Stan cannot handle
%missing cells (for our purposes)
if any(any(ismissing(datatable)))
    error('varargin:missingcells',... %Error code and associated error
    strcat('WARNING: Inputted data have missing cells \n\n',... 
    'There should not be any missing cells in the dataset\n',...
    'Each row in the dataset should contain a measurement\n',...   
    'See ''Preparing Data'' section of UserManual.pdf for more information'));
end

%change verbosity to match true or false
if verbose == 1
    verbosity = false;
elseif verbose == 2
    verbosity = true;
end

%check whether groups, events, or occasions are in the table
if any(strcmpi(colnames,'group'))
    groupnames = unique(datatable.group(:));
    ngroup = length(groupnames);
else
    ngroup = 0;
end

if any(strcmpi(colnames,'event'))
    eventnames = unique(datatable.event(:));
    nevent = length(eventnames);
else
    nevent = 0;
end

if any(strcmpi(colnames,'time'))
    timenames = unique(datatable.time(:));
    ntime = length(timenames);
else 
    ntime = 0;
end


%create a structure array to store information that will later be outputted
REL = struct;
REL.filename = datatable.Properties.Description;
REL.niter = niter;
REL.nchains = nchains;

%determine how cmdstan will be set up
%analysis variable will indicate whether group or events need to be
%considered when sending code to cmdstan
%analysis:
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider
%5 - possible occasions to consider

if (ntime == 0)
    if (ngroup == 0) && (nevent == 0)
        analysis = 1;
    elseif (ngroup > 0) && (nevent == 0)
        analysis = 2;
    elseif (ngroup == 0) && (nevent > 0)
        analysis = 3;
    elseif (ngroup > 0) && (nevent > 0)
        analysis = 4;
    end
elseif (ntime > 0)
    analysis = 5;
end

%show a gui that indicates data are processing in cmdstan if the user
%specified to do so
if showgui == 2
    era_relgui = findobj('Tag','era_relgui');
    if ~isempty(era_relgui)
        close(era_relgui);
    end
    %define parameters for figure position
    figwidth = 400;
    figheight = 150;
    fsize = get(0,'DefaultTextFontSize');

    %define space between rows and first row location
    rowspace = 30;
    row = figheight - rowspace*2;

    %initialize gui
    era_relgui= figure('unit','pix','Visible','off',...
      'position',[400 400 figwidth figheight],...
      'menub','no','numbertitle','off','resize','off');

    movegui(era_relgui,'center');

    %Write text
    uicontrol(era_relgui,'Style','text','fontsize',fsize+6,...
        'HorizontalAlignment','center',...
        'String','Data are processing in CmdStan',...
        'Position',[0 row figwidth 25]);          

    %display gui
    set(era_relgui,'Visible','on');
    
    %tag gui
    era_relgui.Tag = 'era_relgui';
    
    %pause a moment so the gui will be displayed
    pause(.02)

end

switch analysis
    case 1 %no groups or event types to consider
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,'id');
        id2 = zeros(0,length(datatable.id));

        for i = 1:length(datatable.id)
            if i == 1
                id2(1) = 1;
                count = 1;
            elseif i > 1 && strcmp(char(datatable.id(i)), char(datatable.id(i-1)))
                id2(i) = count;
            elseif i > 1 && ~strcmp(datatable.id(i), datatable.id(i-1))
                count = count+1;
                id2(i) = count;
            end
        end

        datatable.id2 = id2(:);
        
        %store the datatable used in REL
        REL.data = datatable;
        REL.analysis = 'ic';
        
        %groups and events were not considered
        REL.groups = 'none';
        REL.events = 'none';
        

        %create cmdstan syntax
        stan_in = {
          'data {' 
          '  int<lower=0> NOBS; //number of obs'
          '  int<lower=0> NSUB; //number of subj'
          '  int<lower=0, upper=NSUB> id[NOBS]; //subject id;'
          '  real meas[NOBS];'
          '}'
          'parameters {'
          '  real mu;'
          '  real<lower=0> sig_u;'
          '  real<lower=0> sig_e;'
          '  vector[NSUB] u_raw;'
          '}'
          'transformed parameters {'
          '  vector[NSUB] u;'
          '  u = mu + sig_u*u_raw;'
          '}'
          'model {'
          '  u_raw ~ normal(0,1);'
          '    meas ~ normal(u[id], sig_e);'
          '  '
          '  mu ~ normal(0,100);'
          '  sig_u ~ cauchy(0,40);'
          '  sig_e ~ cauchy(0,40);'
          '}'
        };
        
        %create data structure to pass to cmdstan
        data = struct(...
            'NOBS',length(datatable.id), ... %number of observations
            'NSUB', length(unique(datatable.id)),... %number of participants
            'id', datatable.id2,... %id variable
            'meas', datatable.meas); %measurement variable
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        modelname = strcat('cmdstan',char(date));
        
        %fit model in cmdstan
        fit = stan('model_code', stan_in,... 
            'model_name', modelname,...
            'data', data,... 
            'iter', niter,...
            'chains', nchains,... 
            'refresh', niter/10,... 
            'verbose', verbosity,... 
            'file_overwrite', true);
        
        %don't let the user continue to use the Matlab command window
        fit.block();
        
        %extract and store cmdstan output
        REL.out.mu = fit.extract('pars','mu').mu;
        REL.out.sig_u = fit.extract('pars','sig_u').sig_u;
        REL.out.sig_e = fit.extract('pars','sig_e').sig_e; 
        
        %label is simply measure (no events or groups to deal with)
        REL.out.labels = 'measure';
        
        REL.out.conv.data = era_storeconv(fit);
        
        
    case 2 %possible multiple groups 

        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'group','id'});
        id2 = zeros(0,length(datatable.id));

        for i = 1:length(datatable.id)
            if i == 1
                id2(1) = 1;
                count = 1;
            elseif i > 1 &&...
                    strcmp(char(datatable.id(i)), char(datatable.id(i-1)))...
                    &&...
                    strcmp(char(datatable.group(i)), char(datatable.group(i-1)))
                id2(i) = count;
            elseif i > 1 &&...
                    ~strcmp(char(datatable.id(i)), char(datatable.id(i-1)))...
                    &&...
                    strcmp(char(datatable.group(i)), char(datatable.group(i-1)))
                count = count+1;
                id2(i) = count;
            elseif i > 1 &&...
                    ~strcmp(char(datatable.id(i)), char(datatable.id(i-1)))...
                    &&...
                    ~strcmp(char(datatable.group(i)), char(datatable.group(i-1)))
                count = 1;
                id2(i) = 1;
            end
        end

        datatable.id2 = id2(:);
        
        %create labels for the groups
        REL.data = datatable;
        REL.analysis = 'ic';
        groupnames = unique(datatable.group(:));
        ngroup = length(groupnames);
        
        %all names are forced to be strings for ease of use later on
        if isnumeric(groupnames)
            groupnames = num2str(groupnames);
        end
        
        %store group names. separate events are not a consideration
        REL.groups = groupnames;
        REL.events = 'none';
        
        
        %create model syntax for cmdstan
        stan_in{1,1} = 'data {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> NG%d; //number of obs in g%d',i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> JG%d; //number of subj in g%d',i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0, upper=JG%d> id_G%d[NG%d]; //subject id g%d;',...
                i,i,i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real meas_G%d[NG%d];',i,i);
        end    

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'parameters {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real mu_G%d;',i);
        end 

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_u_G%d;',i);
        end 

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_e_G%d;',i);
        end 

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  vector[JG%d] u_raw_G%d;',i,i);
        end 

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'transformed parameters {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  vector[JG%d] u_G%d;',i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  u_G%d = mu_G%d + sig_u_G%d*u_raw_G%d;',...
                i,i,i,i);
        end

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'model {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  u_raw_G%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('    meas_G%d ~ normal(u_G%d[id_G%d], sig_e_G%d);',...
                i,i,i,i);
        end

        stan_in{end+1,1} = '  ';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  mu_G%d ~ normal(0,100);',i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  sig_u_G%d ~ cauchy(0,40);',i);
        end
        
        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  sig_e_G%d ~ cauchy(0,40);',i);
        end
        
        stan_in{end+1,1} = '}';
        
        
        %store cmdstan model syntax
        REL.stan_in = stan_in;
        
        %create structure array for the data to be sent to cmdstan
        data = struct;
        
        for i=1:ngroup
            fieldname = sprintf('NG%d',i);
            fieldvalue = height(datatable(ismember(datatable.group,...
                groupnames(i)),1)); 
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:ngroup
            fieldname = sprintf('JG%d',i);
            dummy = datatable(ismember(datatable.group,...
                groupnames(i)),:);
            fieldvalue = length(unique(dummy.id));
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:ngroup
            fieldname = sprintf('id_G%d',i);
            fieldvalue = datatable{ismember(datatable.group,...
                groupnames(i)),{'id2'}};
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:ngroup
            fieldname = sprintf('meas_G%d',i);
            fieldvalue = datatable{ismember(datatable.group,...
                groupnames(i)),{'meas'}};
            data.(fieldname) = fieldvalue;
        end
        
        %execute model code 
        
        modelname = strcat('cmdstan',char(date));
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        %run cmdstan
        fit = stan('model_code', stan_in,... 
            'model_name', modelname,...
            'data', data,...
            'iter', niter,...
            'chains', nchains,...
            'refresh', niter/10,...
            'verbose', verbosity,...
            'file_overwrite', true,...
            'seed', 12345);
        
        %block user from using Matlab command window
        fit.block();
        
        %create fields for storing parsed cmdstan outputs
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};  
        
        %store cmdstand outputs
        for i=1:ngroup
            measname = sprintf('mu_G%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.mu(:,end+1) = measvalue;           
        end
        
        for i=1:ngroup
            measname = sprintf('sig_u_G%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_u(:,end+1) = measvalue;           
        end
        
        for i=1:ngroup
            measname = sprintf('sig_e_G%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_e(:,end+1) = measvalue;           
        end
        
        %store labels
        for i=1:ngroup
            REL.out.labels(:,end+1) = groupnames(i);           
        end
        
        REL.out.conv.data = era_storeconv(fit);
        
    case 3 %possible event types but no groups to consider

        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'id','event'}); 
        
        eventnames = unique(datatable.event(:));
        nevent = length(eventnames);
        
        %ensure all event names are strings (makes life easier later)
        if isnumeric(eventnames)
            eventnames = num2str(eventnames);
        end
        
        %store the event names
        REL.events = eventnames;
        
        %no groups to consider
        REL.groups = 'none';
        REL.analysis = 'ic';
        
        %create data structure for storing measurements for each event
        eventarray = struct;
        eventarray.data = [];
        
        for i = 1:nevent
            
            dummytable = datatable(ismember(datatable.event,...
                char(eventnames(i))),:);
            eventarray.data{end+1} = dummytable;
            
        end
        
        %creating that id2 variable to make cmdstan happy
        for j = 1:nevent
            
            id2 = zeros(0,height(eventarray.data{j}));
            
            for i = 1:height(eventarray.data{j})
                
                if i == 1
                    id2(1) = 1;
                    count = 1;
                elseif i > 1 && strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1)))
                    id2(i) = count;
                elseif i > 1 && ~strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1)))
                    count = count+1;
                    id2(i) = count;
                end
            end
            
            eventarray.data{j}.id2 = id2';
            
        end
        
        %store the data
        REL.data = eventarray.data;

        stan_in{1,1} = 'data {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> NE%d; //number of obs in E%d',i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> JE%d; //number of subj in E%d',i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0, upper=JE%d> id_E%d[NE%d]; //subject id E%d;',...
                i,i,i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real meas_E%d[NE%d];',i,i);
        end    

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'parameters {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real mu_E%d;',i);
        end 

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_u_E%d;',i);
        end 

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_e_E%d;',i);
        end 

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  vector[JE%d] u_raw_E%d;',i,i);
        end 

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'transformed parameters {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  vector[JE%d] u_E%d;',i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  u_E%d = mu_E%d + sig_u_E%d*u_raw_E%d;',...
                i,i,i,i);
        end

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'model {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  u_raw_E%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('    meas_E%d ~ normal(u_E%d[id_E%d], sig_e_E%d);',...
                i,i,i,i);
        end

        stan_in{end+1,1} = '  ';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  mu_E%d ~ normal(0,100);',i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  sig_u_E%d ~ cauchy(0,40);',i);
        end
        
        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  sig_e_E%d ~ cauchy(0,40);',i);
        end
        
        stan_in{end+1,1} = '}';
        
        REL.stan_in = stan_in;
        
        %create structure array for the data to be sent to cmdstan
        data = struct;

        for i=1:nevent %number of observations
            fieldname = sprintf('NE%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = height(dummytable); 
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:nevent %number of participants
            fieldname = sprintf('JE%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = length(unique(dummytable.id));
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:nevent %id variable
            fieldname = sprintf('id_E%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = dummytable.id2;
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:nevent %measurement variable
            fieldname = sprintf('meas_E%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = dummytable.meas;
            data.(fieldname) = fieldvalue;
        end
        
        %execute model code 
        
        modelname = strcat('cmdstan',char(date));
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        %run in cmdstan
        fit = stan('model_code', stan_in,...
            'model_name', modelname,...
            'data', data,... 
            'iter', niter,...
            'chains', nchains,...
            'refresh', niter/10,...
            'verbose', verbosity,...
            'file_overwrite', true,...
            'seed', 12345);
        
        %block user from using Matlab command window 
        fit.block();
    
        
        %create fields to store cmdstan output
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};
        REL.out.conv.data = {};
        
        %parse outputs
        for i=1:nevent
            measname = sprintf('mu_E%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.mu(:,end+1) = measvalue;           
        end
        
        for i=1:nevent
            measname = sprintf('sig_u_E%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_u(:,end+1) = measvalue;           
        end
        
        for i=1:nevent
            measname = sprintf('sig_e_E%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_e(:,end+1) = measvalue;           
        end
       
        %labels are the event names
        for i=1:nevent
            REL.out.labels(:,end+1) = eventnames(i);           
        end
       
        REL.out.conv.data = era_storeconv(fit);
        
    case 4 %possible groups and event types to consider

        %to decrease memory load, cmdstan will be run separately for each
        %event. Groups will be processed in the same cmdstan run.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'group','id','event'});
        eventnames = unique(datatable.event(:));
        if isnumeric(eventnames)
            eventnames = num2str(eventnames);
        end
        
        %store event names
        nevent = length(eventnames);
        REL.events = eventnames;
       
        groupnames = unique(datatable.group(:));
        if isnumeric(groupnames)
            groupnames = num2str(groupnames);
        end
        
        %store group names
        ngroup = length(groupnames);
        REL.groups = groupnames;
        
        %create structure for storing data
        eventarray = struct;
        eventarray.data = [];
        
        %store the data for each event separately
        for i = 1:nevent
            
            dummytable = datatable(ismember(datatable.event,...
                char(eventnames(i))),:);
            eventarray.data{end+1} = dummytable;
            
        end
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement. This is
        %done separately for each event since the event data will be passed
        %separately to cmdstan
        
        for j = 1:nevent
            
            id2 = zeros(0,height(eventarray.data{j}));
            
            for i = 1:height(eventarray.data{j})
                if i == 1
                    id2(1) = 1;
                    count = 1;
                elseif i > 1 &&...
                        strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1))) &&...
                        strcmp(char(eventarray.data{j}.group(i)),...
                        char(eventarray.data{j}.group(i-1)))
                    id2(i) = count;
                elseif i > 1 &&...
                        ~strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1))) &&...
                        strcmp(char(eventarray.data{j}.group(i)),...
                        char(eventarray.data{j}.group(i-1)))
                    count = count+1;
                    id2(i) = count;
                elseif i > 1 &&...
                        ~strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1))) &&...
                        ~strcmp(char(eventarray.data{j}.group(i)),...
                        char(eventarray.data{j}.group(i-1)))
                    count = 1;
                    id2(i) = 1;
                end
            end
            
            eventarray.data{j}.id2 = id2';
            
        end
        
        %put the compiled data into REL structure
        REL.data = eventarray.data;
        REL.analysis = 'ic';
        
        %create the fields for storing the parsed cmdstan outputs 
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};
        REL.out.conv.data = {};
        REL.stan_in = {};
        
        %separate syntax will need to be generated for each event because
        %the number of samples will be different
        
        %create syntax for cmdstan
        for j=1:nevent
            
            %create string to be printed and potentially viewed in gui
            str = ['Working on event ' num2str(j) ' of ' num2str(nevent)];

            %print to screen to notify user of progress
            fprintf(strcat('\n\n',str,'\n\n'));
            
            if showgui == 2
                
                %find whether gui exists (the user may have closed it)
                era_relgui = findobj('Tag','era_relgui');

                if ~isempty(era_relgui)

                    %Write text
                    uicontrol(era_relgui,'Style','text',...
                        'fontsize',fsize+6,...
                        'HorizontalAlignment','center',...
                        'String',str,...
                        'Position',[0 row-rowspace*1.5 figwidth 25]);
                    
                    %pause a moment so the gui will be displayed
                    pause(.02)
                end
            end
                     
            clear stan_in
            stan_in{1,1} = 'data {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  int<lower=0> NG%d; //number of obs in g%d',i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  int<lower=0> JG%d; //number of subj in g%d',i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  int<lower=0, upper=JG%d> id_G%d[NG%d]; //subject id g%d;',...
                    i,i,i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real meas_G%d[NG%d];',i,i);
            end    

            stan_in{end+1,1} = '}';
            stan_in{end+1,1} = 'parameters {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real mu_G%d;',i);
            end 

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real<lower=0> sig_u_G%d;',i);
            end 

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real<lower=0> sig_e_G%d;',i);
            end 

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  vector[JG%d] u_raw_G%d;',i,i);
            end 

            stan_in{end+1,1} = '}';
            stan_in{end+1,1} = 'transformed parameters {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  vector[JG%d] u_G%d;',i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  u_G%d = mu_G%d + sig_u_G%d*u_raw_G%d;',...
                    i,i,i,i);
            end

            stan_in{end+1,1} = '}';
            stan_in{end+1,1} = 'model {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  u_raw_G%d ~ normal(0,1);',i);
                stan_in{end+1,1} = ...
                    sprintf('    meas_G%d ~ normal(u_G%d[id_G%d], sig_e_G%d);',...
                    i,i,i,i);
            end

            stan_in{end+1,1} = '  ';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  mu_G%d ~ normal(0,100);',i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  sig_u_G%d ~ cauchy(0,40);',i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  sig_e_G%d ~ cauchy(0,40);',i);
            end

            stan_in{end+1,1} = '}';
            
            %store the cmdstan syntax
            REL.stan_in{end+1} = stan_in;

            %create structure for the data to be sent to cmdstan
            data = struct;

            %prep data structure for cmdstan
            for i=1:ngroup
                fieldname = sprintf('NG%d',i);
                fieldvalue = length(eventarray.data{j}.group(ismember(eventarray.data{j}.group,...
                    groupnames(i)),1)); 
                data.(fieldname) = fieldvalue;
            end

            for i=1:ngroup
                fieldname = sprintf('JG%d',i);
                dummy = eventarray.data{j}(ismember(eventarray.data{j}.group,...
                    groupnames(i)),:);
                fieldvalue = length(unique(dummy.id));
                data.(fieldname) = fieldvalue;
            end

            for i=1:ngroup
                fieldname = sprintf('id_G%d',i);
                fieldvalue = eventarray.data{j}{ismember(eventarray.data{j}.group,...
                    groupnames(i)),{'id2'}};
                data.(fieldname) = fieldvalue;
            end

            for i=1:ngroup
                fieldname = sprintf('meas_G%d',i);
                fieldvalue = eventarray.data{j}{ismember(eventarray.data{j}.group,...
                    groupnames(i)),{'meas'}};
                data.(fieldname) = fieldvalue;
            end

            %execute model code 

            modelname = strcat('cmdstan',char(date));

            fprintf('\nModel is being run in cmdstan\n');
            fprintf('\nThis may take a while depending on the amount of data\n');
            
            %run cmdstan
            fit = stan('model_code', stan_in,... 
                'model_name', modelname,...
                'data', data,... 
                'iter', niter,...
                'chains', nchains,...
                'refresh', niter/10,...
                'verbose', verbosity,...
                'file_overwrite', true,...
                'seed', 12345);
            
            %block user from using Matlab command window
            fit.block();
            
            %parse cmdstan outputs
            for i=1:ngroup
                measname = sprintf('mu_G%d',i);
                measvalue = fit.extract('pars',measname).(measname);
                REL.out.mu(:,end+1) = measvalue;           
            end

            for i=1:ngroup
                measname = sprintf('sig_u_G%d',i);
                measvalue = fit.extract('pars',measname).(measname);
                REL.out.sig_u(:,end+1) = measvalue;           
            end

            for i=1:ngroup
                measname = sprintf('sig_e_G%d',i);
                measvalue = fit.extract('pars',measname).(measname);
                REL.out.sig_e(:,end+1) = measvalue;           
            end

            %store labels for each column by separating the event name and
            %group name with an underscore
            for i=1:ngroup
                REL.out.labels(:,end+1) = strcat(eventnames(j),...
                    '_;_',groupnames(i));           
            end

            REL.out.conv.data{end+1} = era_storeconv(fit);
            
        end %for j=1:nevent
        
    case 5
        
        fprintf('\nPreparing data for analysis...\n');
        
        if (ngroup == 0) && (nevent == 0)
            datatable = sortrows(datatable,{'id','time'});
        elseif (ngroup > 0) && (nevent == 0)
            datatable = sortrows(datatable,{'group','id','time'});
        elseif (ngroup == 0) && (nevent > 0)
            datatable = sortrows(datatable,{'id','time','event'});
        elseif (ngroup > 0) && (nevent > 0)
            datatable = sortrows(datatable,{'group','id','time','event'});
        end
        
        if nevent > 0
            REL.events = eventnames;
        elseif nevent == 0
            REL.events = 'none';
        end
        
        if ngroup > 0
            REL.groups = groupnames;
        elseif ngroup == 0
            REL.groups = 'none';
        end
       
        REL.time = timenames;
  
        %parse data based on number of events and groups
        %create structure for storing data
        darray = struct;
        darray.data = [];
        darray.names = [];
        
        if (ngroup == 0) && (nevent == 0)
            darray.data{end+1} = datatable; 
            darray.names{end+1} = 'none';
        elseif (ngroup > 0) && (nevent == 0)
            for i = 1:ngroup
                dummytable = datatable(ismember(datatable.group,...
                    char(groupnames(i))),:);
                darray.data{end+1} = dummytable;
                darray.names{end+1} = char(groupnames(i));
            end
        elseif (ngroup == 0) && (nevent > 0)
            for i = 1:nevent
                dummytable = datatable(ismember(datatable.event,...
                    char(eventnames(i))),:);
                darray.data{end+1} = dummytable;
                darray.names{end+1} = char(eventnames(i));
            end
        elseif (ngroup > 0) && (nevent > 0)
            for i = 1:ngroup
                dummytable = datatable(ismember(datatable.group,...
                    char(groupnames(i))),:);
                for j = 1:nevent
                    dummytable2 = dummytable(ismember(dummytable.event,...
                        char(eventnames(j))),:);
                    darray.data{end+1} = dummytable2;
                    darray.names{end+1} = strcat(char(groupnames(i)),...
                        '_;_',char(eventnames(j)));
                end
            end
        end
        
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement. 
        
        for da = 1:length(darray.names)
            warray = darray.data{da};
            id2 = zeros(0,height(warray));
            time2 = zeros(0,height(warray));
            trl2 = zeros(0,height(warray));
            for i = 1:height(warray)
                %recode ids
                if i == 1
                    id2(1) = 1;
                    count = 1;
                elseif i > 1 && strcmp(char(warray.id(i)),...
                        char(warray.id(i-1)))
                    id2(i) = count;
                elseif i > 1 && ~strcmp(char(warray.id(i)),...
                        char(warray.id(i-1)))
                    count = count+1;
                    id2(i) = count;
                end
                
                %recode time
                time2(i) = find(strcmp(timenames,warray.time(i)));
                
                %recode trials
                if i == 1
                    trl2(1) = 1;
                    trlcount = 2;
                elseif i > 1 && strcmp(char(warray.id(i)),...
                        char(warray.id(i-1))) && ...
                        strcmp(char(warray.time(i)),...
                        char(warray.time(i-1)))
                    trl2(i) = trlcount;
                    trlcount = trlcount + 1;
                elseif i > 1 && strcmp(char(warray.id(i)),...
                        char(warray.id(i-1))) && ...
                        ~strcmp(char(warray.time(i)),...
                        char(warray.time(i-1)))
                    trl2(i) = 1;
                    trlcount = 2;
                elseif i > 1 && ~strcmp(char(warray.id(i)),...
                        char(warray.id(i-1)))
                    trl2(i) = 1;
                    trlcount = 2;
                end
            end
            warray.id2 = id2';
            warray.time2 = time2';
            warray.trl2 = trl2';
            darray.data{da} = warray;
        end
        
        for da = 1:length(darray.names)
            warray = darray.data{da};
             for i = 1:height(warray)
                if i == 1
                    idxtrl_count = 1;
                    idxtrl_count_base = 1;
                    
                    idxtrl = idxtrl_count;
                    poss_idxtrl = idxtrl_count;
                    idxtrl_count = idxtrl_count + 1;
                elseif i > 1 
                    if warray.id2(i) == warray.id2(i-1)
                        if warray.trl2(i) == (warray.trl2(i-1)+1)
                            poss_idxtrl(end+1) = idxtrl_count;
                            idxtrl(end+1) = idxtrl_count;                            
                            idxtrl_count = idxtrl_count + 1;
                        elseif warray.trl2(i) ~= (warray.trl2(i-1)+1)
                            idxtrl_count = idxtrl_count_base;
                            poss_idxtrl(end+1) = idxtrl_count;
                            idxtrl(end+1) = idxtrl_count;
                            idxtrl_count = idxtrl_count + 1;
                        end
                    elseif warray.id2(i) ~= warray.id2(i-1)
                        idxtrl_count_base = max(poss_idxtrl) + 1;
                        idxtrl_count = idxtrl_count_base;
                        poss_idxtrl = idxtrl_count;
                        idxtrl(end+1) = idxtrl_count;
                        idxtrl_count = idxtrl_count + 1;
                    end
                end
            end
            
            for i = 1:height(warray)
                if i == 1
                    idxtim_count = 1;
                    idxtim = idxtim_count;
                elseif i > 1 
                    if warray.id2(i) == warray.id2(i-1)
                        if warray.time2(i) == (warray.time2(i-1))
                            idxtim(end+1) = idxtim_count;    
                        elseif warray.time2(i) ~= (warray.time2(i-1))
                            idxtim_count = idxtim_count + 1;
                            idxtim(end+1) = idxtim_count;
                        end
                    elseif warray.id2(i) ~= warray.id2(i-1)
                        idxtim_count = idxtim_count + 1;
                        idxtim(end+1) = idxtim_count;
                    end
                end
            end
            
            for i = 1:height(warray)
                if i == 1
                    trlxtim_count = 1;
                    ind_trlxtim = [unique(warray.time2) ...
                        zeros(1,length(unique(warray.time2)))'];
                    
                    for j = 1:length(unique(warray.time2))
                        if j == 1
                            ind_trlxtim(j,2) = 1;
                            maxprev = max(warray{warray.time2 == 1,'trl2'});
                        elseif j > 1
                            ind_trlxtim(j,2) = maxprev + 1;
                            maxprev = maxprev + 1 +... 
                                max(warray{warray.time2 == 1,'trl2'});
                        end                      
                    end
                    
                    trlxtim = trlxtim_count;
                    trlxtim_count = trlxtim_count + 1;
                elseif i > 1 
                    if warray.time2(i) == warray.time2(i-1)
                        if warray.trl2(i) == (warray.trl2(i-1)+1)                           
                            trlxtim(end+1) = trlxtim_count;                            
                            trlxtim_count = trlxtim_count + 1;
                        elseif warray.trl2(i) ~= (warray.trl2(i-1)+1)
                            trlxtim_count = ind_trlxtim(warray.time2(i),2);
                            trlxtim(end+1) = trlxtim_count;
                            trlxtim_count = trlxtim_count + 1;
                        end
                    elseif warray.time2(i) ~= warray.time2(i-1)
                        trlxtim_count = ind_trlxtim(warray.time2(i),2);
                        trlxtim(end+1) = trlxtim_count;
                        trlxtim_count = trlxtim_count + 1;
                    end
                end
            end
            warray.idxtrl = idxtrl';
            warray.idxtim = idxtim';
            warray.trlxtim = trlxtim';
            darray.data{da} = warray;
        end
        
        %put the compiled data into REL structure
        REL.data = datatable;
        REL.analysis = 'trt';
        
        if showgui == 2
            
            %find whether gui exists (the user may have closed it)
            era_relgui = findobj('Tag','era_relgui');
            
            if ~isempty(era_relgui)
                
                %Write text
                uicontrol(era_relgui,'Style','text',...
                    'fontsize',fsize+6,...
                    'HorizontalAlignment','center',...
                    'String','Data are crunching. This will take a long time!',...
                    'Position',[0 row-rowspace*1.5 figwidth 25]);
                
                %pause a moment so the gui will be displayed
                pause(.02)
            end
        end
        
        %define number of data chunks to crunch through
        ndchunks = length(darray.names);
        
        REL.stan_in = [];
        
        for i = 1:ndchunks
            
            %create string to be printed and potentially viewed in gui
            str = ['Working on event/group ' num2str(i) ' of ' num2str(ndchunks)];
            
            %print to screen to notify user of progress
            fprintf(strcat('\n\n',str,'\n\n'));
            
            if showgui == 2
                
                %find whether gui exists (the user may have closed it)
                era_relgui = findobj('Tag','era_relgui');
                
                if ~isempty(era_relgui)
                    
                    %Write text
                    uicontrol(era_relgui,'Style','text',...
                        'fontsize',fsize+6,...
                        'HorizontalAlignment','center',...
                        'String',str,...
                        'Position',[0 row-rowspace*1.5 figwidth 25]);
                    
                    %pause a moment so the gui will be displayed
                    pause(.02)
                end
            end
            
            stan_in = {
                'data { '
                '  int<lower=1> NOBS; //total number of observations '
                '  int<lower=1> NSUB; //total number of subjects'
                '  int<lower=1> NOCC; //total number of occasions'
                '  int<lower=1> NTRL; //total number of trials'
                '  int<lower=1> NTID; // total unique for trial*id'
                '  int<lower=1> NOID; // total unique for occ*id'
                '  int<lower=1> NTO; // total unique for trial*occ'
                '  int<lower=1, upper=NSUB> id[NOBS]; //id variable'
                '  int<lower=1, upper=NOCC> occ[NOBS]; //occ variable'
                '  int<lower=1, upper=NTRL> trl[NOBS]; //trl variable'
                '  int<lower=1, upper=NTID> trlxid[NOBS]; //trlxid variable'
                '  int<lower=1, upper=NOID> occxid[NOBS]; //occxid variable'
                '  int<lower=1, upper=NTO> trlxocc[NOBS]; //trlxocc variable'
                '  vector[NOBS] meas; //measurements'
                '}'
                'parameters {'
                '  real pop_int; //population intercept'
                '  real<lower=0> sig_id; //id-level std dev'
                '  real<lower=0> sig_occ; //occ-level std dev'
                '  real<lower=0> sig_trl; //trl-level std dev'
                '  real<lower=0> sig_err; //residual std dev'
                '  real<lower=0> sig_trlxid; //trlxid std dev'
                '  real<lower=0> sig_occxid; //occxid std dev'
                '  real<lower=0> sig_trlxocc; //trlxocc std dev'
                '  vector[NSUB] id_raw; //id means'
                '  vector[NOCC] occ_raw; //occ means'
                '  vector[NTRL] trl_raw; //trl means'
                '  vector[NTID] trlxid_raw; //trlxid means'
                '  vector[NOID] occxid_raw; //occxid means'
                '  vector[NTO] trlxocc_raw; //trlxocc means'
                '}'
                'transformed parameters {'
                '  vector[NSUB] id_terms; //id-level terms'
                '  vector[NOCC] occ_terms; //occ-level terms'
                '  vector[NTRL] trl_terms; //trl-level terms'
                '  vector[NTID] trlxid_terms; //trlxid terms'
                '  vector[NOID] occxid_terms; //trlxocc terms'
                '  vector[NTO] trlxocc_terms; //trlxocc terms'
                '  id_terms = sig_id * id_raw;'
                '  occ_terms = sig_occ * occ_raw;'
                '  trl_terms = sig_trl * trl_raw;'
                '  trlxid_terms = sig_trlxid * trlxid_raw;'
                '  occxid_terms = sig_occxid * occxid_raw;'
                '  trlxocc_terms = sig_trlxocc * trlxocc_raw;'
                '}'
                'model {'
                '  vector[NOBS] y_hat;'
                '    y_hat = pop_int + id_terms[id] + '
                '    occ_terms[occ] + trl_terms[trl] +'
                '    trlxid_terms[trlxid] + occxid_terms[occxid] + '
                '    trlxocc_terms[trlxocc];'
                '  meas ~ normal(y_hat,sig_err);'
                '  id_raw ~ normal(0,1);'
                '  occ_raw ~ normal(0,1);'
                '  trl_raw ~ normal(0,1);'
                '  trlxid_raw ~ normal(0,1);'
                '  occxid_raw ~ normal(0,1);'
                '  trlxocc_raw ~ normal(0,1);'
                '  sig_id ~ cauchy(0,10);'
                '  sig_occ ~ cauchy(0,.05);'
                '  sig_trl ~ cauchy(0,10);'
                '  sig_err ~ cauchy(0,20);'
                '  sig_trlxid ~ cauchy(0,5);'
                '  sig_occxid ~ cauchy(0,5);'
                '  sig_trlxocc ~ cauchy(0,.05);'
                '}'
                };
            
            
            %store the cmdstan syntax
            REL.stan_in{end+1} = stan_in;
            
            %create structure for the data to be sent to cmdstan
            data = struct;
            
            %prep data structure for cmdstan
            
            fieldname = 'NOBS';
            fieldvalue = height(darray.data{i});
            data.(fieldname) = fieldvalue;
            
            fieldname = 'NOCC';
            fieldvalue = length(unique(darray.data{i}.time));
            data.(fieldname) = fieldvalue;
            
            
            fieldname = 'NSUB';
            fieldvalue = length(unique(darray.data{i}.id2));
            data.(fieldname) = fieldvalue;
            
            
            fieldname = 'NTRL';
            fieldvalue = length(unique(darray.data{i}.trl2));
            data.(fieldname) = fieldvalue;
            
            fieldname = 'occ';
            fieldvalue = darray.data{i}.time2;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'id';
            fieldvalue = darray.data{i}.id2;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'trl';
            fieldvalue = darray.data{i}.trl2;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'meas';
            fieldvalue = darray.data{i}.meas;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'trlxid';
            fieldvalue = darray.data{i}.idxtrl;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'occxid';
            fieldvalue = darray.data{i}.idxtim;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'trlxocc';
            fieldvalue = darray.data{i}.trlxtim;
            data.(fieldname) = fieldvalue;
            
            fieldname = 'NTID';
            fieldvalue = length(unique(darray.data{i}.idxtrl));
            data.(fieldname) = fieldvalue;
            
            fieldname = 'NOID';
            fieldvalue = length(unique(darray.data{i}.idxtim));
            data.(fieldname) = fieldvalue;
            
            fieldname = 'NTO';
            fieldvalue = length(unique(darray.data{i}.trlxtim));
            data.(fieldname) = fieldvalue;
            
            %execute model code
            
            modelname = strcat('cmdstan',char(date));
            
            fprintf('\nModel is being run in cmdstan\n');
            fprintf('\nThis may take a while depending on the amount of data\n');
            
            stan_control = struct;
            stan_control.delta = .98;
            stan_control.t0 = 15;
            
            %run cmdstan
            fit = stan('model_code', stan_in,...
                'model_name', modelname,...
                'data', data,...
                'iter', niter,...
                'chains', nchains,...
                'refresh', niter/20,...
                'verbose', verbosity,...
                'file_overwrite', true,...
                'control', stan_control,...
                'seed', 12345);
            
            %block user from using Matlab command window
            fit.block();
            
            %what happens if niter/20 is not an evenly divisly number of niter
            %or decimal number?
            
            if i == 1
                %create the fields for storing the parsed cmdstan outputs
                REL.out = [];
                REL.out.mu = [];
                REL.out.sig_id = [];
                REL.out.sig_occ = [];
                REL.out.sig_trl = [];
                REL.out.sig_trlxid = [];
                REL.out.sig_occxid = [];
                REL.out.sig_trlxocc =[];
                REL.out.sig_err = [];
                REL.out.labels = {};
                REL.out.conv.data = {};
            end
            
            measname ='pop_int';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.mu(:,end+1) = measvalue;
            
            measname = 'sig_id';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_id(:,end+1) = measvalue;
            
            measname = 'sig_occ';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_occ(:,end+1) = measvalue;
            
            measname = 'sig_trl';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_trl(:,end+1) = measvalue;
            
            measname = 'sig_trlxid';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_trlxid(:,end+1) = measvalue;
            
            measname = 'sig_occxid';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_occxid(:,end+1) = measvalue;
            
            measname = 'sig_trlxocc';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_trlxocc(:,end+1) = measvalue;
            
            measname = 'sig_err';
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_err(:,end+1) = measvalue;
            
            REL.out.conv.data{end+1} = era_storeconv(fit,1);
            
        end
        
        REL.out.labels = darray.names;
        
end %switch analysis

%close the gui if one was shown
if showgui == 2
    era_relgui = findobj('Tag','era_relgui');
    if ~isempty(era_relgui)
        close(era_relgui);
    end
end

REL.eraver = eraver;
RELout = REL;

end

function convstats = era_storeconv(fit,trt)
%pull out r_hats from Stanfit model

if nargin < 2
    trt = 0;
end

%pull summary using the print function
%there's no other way to get at the r_hat values, so output will also be
%printed in the command window
output = print(fit);

%define cell array that holds the r_hats
convstats = cell(0,3);
convstats{1,1} = 'name';
convstats{1,2} = 'n_eff';
convstats{1,3} = 'r_hat';

if trt == 0
    %cycle through the important inputs
    %depending on whether there are multiple events/groups there may be
    %multiple mu_, sig_u_, and sig_e_
    inp2check = {'lp__' 'mu_' 'sig_u_' 'sig_e_'};
    
    for j = 1:length(inp2check)
        check = strncmp(output,inp2check{j},length(inp2check{j}));
        ind2check = find(check == 1);
        for i = 1:length(ind2check)
            pullrow = strsplit(output{ind2check(i),:},' ');
            label = pullrow{1};
            neff = str2num(pullrow{end-2});
            rhat = str2num(pullrow{end});
            convstats{end+1,1} = label;
            convstats{end,2} = neff;
            convstats{end,3} = rhat;
        end
    end
elseif trt == 1
    inp2check = {'lp__'};
    
    inp2check{end+1} = 'pop_int';
    inp2check{end+1} = 'sig_id';
    inp2check{end+1} = 'sig_occ';
    inp2check{end+1} = 'sig_trl';
    inp2check{end+1} = 'sig_trlxid';
    inp2check{end+1} = 'sig_occxid';
    inp2check{end+1} = 'sig_trlxocc';
    inp2check{end+1} = 'sig_err';
    
    for j = 1:length(inp2check)
        check = strncmp(output,inp2check{j},length(inp2check{j}));
        ind2check = find(check == 1);
        for i = 1:length(ind2check)
            pullrow = strsplit(output{ind2check(i),:},' ');
            label = pullrow{1};
            neff = str2num(pullrow{end-2});
            rhat = str2num(pullrow{end});
            convstats{end+1,1} = label;
            convstats{end,2} = neff;
            convstats{end,3} = rhat;
        end
    end
end
end
