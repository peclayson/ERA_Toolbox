function RELout = era_computerel(varargin)
%Prepare and execute cmdstan code for dependability analyses
%
%era_computerel('data',era_datatable,'chains',3,'iter',1000)
%
%Lasted Updated 2/25/19
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
%2/25/19 PC
% added functionality for computing test-retest reliability

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
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', verbosity, 'file_overwrite', true);
        
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
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', verbosity, 'file_overwrite', true);
        
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
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', verbosity, 'file_overwrite', true);
        
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
            fit = stan('model_code', stan_in, 'model_name', modelname,...
                'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
                niter/10, 'verbose', verbosity, 'file_overwrite', true);
            
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

        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement. 
        
        id2 = zeros(0,height(datatable));
        time2 = zeros(0,height(datatable));
        trl2 = zeros(0,height(datatable));
        
        for i = 1:height(datatable)
            %recode ids
            if i == 1
                id2(1) = 1;
                count = 1;
            elseif i > 1 && strcmp(char(datatable.id(i)),...
                    char(datatable.id(i-1)))
                id2(i) = count;
            elseif i > 1 && ~strcmp(char(datatable.id(i)),...
                    char(datatable.id(i-1)))
                count = count+1;
                id2(i) = count;
            end
            
            %recode time
            time2(i) = find(strcmp(timenames,datatable.time(i)));
            
            if nevent == 0
                %recode trials
                if i == 1
                    trl2(1) = 1;
                    trlcount = 2;
                elseif i > 1 && strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1))) && ... 
                        strcmp(char(datatable.time(i)),...
                        char(datatable.time(i-1)))
                    trl2(i) = trlcount;
                    trlcount = trlcount + 1;
                elseif i > 1 && strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1))) && ... 
                        ~strcmp(char(datatable.time(i)),...
                        char(datatable.time(i-1)))
                    trl2(i) = 1;
                    trlcount = 2;
                elseif i > 1 && ~strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1)))
                    trl2(i) = 1;
                    trlcount = 2;
                end
            elseif nevent > 0
                %recode trials
                if i == 1
                    trl2(1) = 1;
                    trlcount = 2;
                elseif i > 1 && strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1))) && ... 
                        strcmp(char(datatable.event(i)),...
                        char(datatable.event(i-1))) && ... 
                        strcmp(char(datatable.time(i)),...
                        char(datatable.time(i-1)))
                    trl2(i) = trlcount;
                    trlcount = trlcount + 1;
                elseif (i > 1 && strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1))) && ... 
                        ~strcmp(char(datatable.event(i)),...
                        char(datatable.event(i-1)))) || ...
                        (i > 1 && strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1))) && ... 
                        strcmp(char(datatable.event(i)),...
                        char(datatable.event(i-1))) && ...
                        strcmp(char(datatable.time(i)),...
                        char(datatable.time(i-1))))
                    trl2(i) = 1;
                    trlcount = 2;
                elseif i > 1 && ~strcmp(char(datatable.id(i)),...
                        char(datatable.id(i-1))) 
                    trl2(i) = 1;
                    trlcount = 2;
                end
            end
        end
        
        datatable.id2 = id2';
        datatable.time2 = time2';
        datatable.trl2 = trl2';
        
        
        
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
                    dummytable2 = datatable(ismember(dummytable.event,...
                        char(eventnames(j))),:);
                    darray.data{end+1} = dummytable2;
                    darray.names{end+1} = strcat(char(groupnames(j)),...
                        '_',char(eventnames(i)));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%check the
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%indexing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%here.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%need to
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%do this
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%of the
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%interaction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%variables
        idxtrl_count = 1;
        idxtim = 1;
        trlxtim = 1;
        
        for da = 1:length(darray.names)
            warray = darray.data{da};
            for i = 1:height(warray)
                if i == 1
                    idxtrl = idxtrl_count;
                    idxtrl_count = idxtrl_count + 1;
                    idxtrl_count_base = 1;
                    poss_idxtrl = 1;
                elseif i > 1 
                    if warray.id2(i) == warray.id2(i-1)
                        if warray.trl2(i) == (warray.trl2(i-1)+1)
                            idxtrl(end+1) = idxtrl_count;
                            poss_idxtrl(end+1) = idxtrl_count;
                            idxtrl_count = idxtrl_count + 1;
                        elseif warray.trl2(i) ~= (warray.trl2(i-1)+1)
                            poss_idxtrl(end+1) = idxtrl_count;
                            idxtrl_count = idxtrl_count_base;
                            idxtrl(end+1) = idxtrl_count;
                            idxtrl_count = idxtrl_count + 1;
                        end
                    elseif warray.id2(i) ~= warray.id2(i-1)
                        idxtrl_count_base = max(poss_idxtrl) + 1;
                        idxtrl_count = idxtrl_count_base;
                        poss_idxtrl = [];
                        idxtrl(end+1) = idxtrl_count;
                        idxtrl_count = idxtrl_count + 1;
                    end
                end
            end
            warray.idxtrl = idxtrl';
            darray.data{da} = warray;
        end
        
        %put the compiled data into REL structure
        REL.data = datatable;
        
        %create the fields for storing the parsed cmdstan outputs 
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_id = [];
        REL.out.sig_tim = [];
        REL.out.sig_trl = [];
        REL.out.sig_idxtrl = [];
        REL.out.sig_idxtim = [];
        REL.out.sig_trlxtim =[];
        REL.out.sig_e = [];
        REL.out.labels = {};
        REL.out.conv.data = {};
        REL.stan_in = {};
        
        %parse data into chunks based on event and group
        
        
        
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
        
        clear stan_in
        stan_in{1,1} = 'data {';
        
        for i = 1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NOBS%d; //number of obs in chunk %d',i,i); %#ok<*AGROW>
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NSUB%d; //number of subj in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NOCC%d; //number of occ in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NTRL%d; //number of trials in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NTID%d; //number of unique trlxid in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NOID%d; //number of unique occxid in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NTO%d; //number of unique trlxocc in chunk %d',i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i = 1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NSUB%d> id%d[NOBS%d]; //id variable in chunk %d',...
                i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NOCC%d> occ%d[NOBS%d]; //occ variable in chunk %d',...
                i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NTRL%d> trl%d[NOBS%d]; //trl variable in chunk %d',...
                i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NTID%d> trlxid%d[NOBS%d]; //trlxid variable in chunk %d',...
                i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NOID%d> occxid%d[NOBS%d]; //occxid variable in chunk %d',...
                i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NTO%d> trlxocc%d[NOBS%d]; //trlxocc variable in chunk %d',...
                i,i,i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i = 1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NOBS%d] meas%d; //measurements in chunk %d',i,i,i);
        end
        
        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = '  ';
        stan_in{end+1,1} = 'parameters {';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  real pop_int%d; //population intercept in chunk %d',i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_id%d; //id-level std dev in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_occ%d; //occ-level std dev in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_trl%d; //trl-level std dev in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_err%d; //residual std dev in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_trlxid%d; //trlxid std dev in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_occxid%d; //occxid std dev in chunk %d',i,i);
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_trlxocc%d; //trlxocc std dev in chunk %d',i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NSUB%d] id_raw%d; //id means in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NOCC%d] occ_raw%d; //occ means in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NTRL%d] trl_raw%d; //trl means in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NTID%d] trlxid_raw%d; //trlxid means in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NOID%d] occxid_raw%d; //occxid means in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NTO%d] trlxocc_raw%d; //trlxocc means in chunk %d',i,i,i);
        end
        
        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = '  ';
        stan_in{end+1,1} = 'transformed parameters {';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NSUB%d] id_terms%d; //id-level terms in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NOCC%d] occ_terms%d; //occ-level terms in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NTRL%d] trl_terms%d; //trl-level terms in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NTID%d] trlxid_terms%d; //trlxid terms in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NOID%d] occxid_terms%d; //occxid terms in chunk %d',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  vector[NTO%d] trlxocc_terms%d; //trlxocc terms in chunk %d',i,i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  id_terms%d = sig_id%d * id_raw%d;',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  occ_terms%d = sig_occ%d * occ_raw%d;',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  trl_terms%d = sig_trl%d * trl_raw%d;',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  trlxid_terms%d = sig_trlxid%d * trlxid_raw%d;',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  occxid_terms%d = sig_occxid%d * occxid_raw%d;',i,i,i);
            stan_in{end+1,1} = ...
                sprintf('  trlxocc_terms%d = sig_trlxocc%d * trlxocc_raw%d;',i,i,i);
        end
        
        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = '  ';
        stan_in{end+1,1} = 'model {';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NOBS%d] y_hat%d;',i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  for(i in 1:NOBS%d)',i);
            stan_in{end+1,1} = ...
                sprintf('    y_hat%d[i] = pop_int%d + id_terms%d[id%d[i]] +',i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('    occ_terms%d[occ%d[i]] + trl_terms%d[trl%d[i]] +',i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('    trlxid_terms%d[trlxid%d[i]] + occxid_terms%d[occxid%d[i]] +',i,i,i,i);
            stan_in{end+1,1} = ...
                sprintf('    trlxocc_terms%d[trlxocc%d[i]];',i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  meas%d ~ normal(y_hat%d,sig_err%d);',i,i,i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  id_raw%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  occ_raw%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  trl_raw%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  trlxid_raw%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  occxid_raw%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  trlxocc_raw%d ~ normal(0,1);',i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  sig_id%d ~ cauchy(0,20);',i);
            stan_in{end+1,1} = ...
                sprintf('  sig_occ%d ~ cauchy(0,.1);',i);
            stan_in{end+1,1} = ...
                sprintf('  sig_trl%d ~ cauchy(0,20);',i);
            stan_in{end+1,1} = ...
                sprintf('  sig_err%d ~ cauchy(0,20);',i);
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  sig_trlxid%d ~ cauchy(0,10);',i);
            stan_in{end+1,1} = ...
                sprintf('  sig_occxid%d ~ cauchy(0,10);',i);
            stan_in{end+1,1} = ...
                sprintf('  sig_trlxocc%d ~ cauchy(0,.05);',i);
        end
        
        stan_in{end+1,1} = '}';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
        %define number of data chunks to crunch through
        ndchunks = length(darray.names);
        clear stan_in
        stan_in{1,1} = 'data {';
        
        for i = 1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NOBS;'); %#ok<*AGROW>
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NSUB;');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NOCC;');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NTRL;');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NTID;');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NOID;');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1> NTO;');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i = 1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NSUB> id[NOBS];');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NOCC> occ[NOBS];');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NTRL> trl[NOBS];');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NTID> trlxid[NOBS];');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NOID> occxid[NOBS];');
            stan_in{end+1,1} = ...
                sprintf('  int<lower=1, upper=NTO> trlxocc[NOBS];');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i = 1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NOBS] meas;');
        end
        
        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = '  ';
        stan_in{end+1,1} = 'parameters {';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  real pop_int;');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_id;');
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_occ;');
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_trl;');
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_err;');
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_trlxid;');
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_occxid;');
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_trlxocc;');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NSUB] id_raw;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NOCC] occ_raw;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NTRL] trl_raw;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NTID] trlxid_raw;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NOID] occxid_raw;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NTO] trlxocc_raw;');
        end
        
        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = '  ';
        stan_in{end+1,1} = 'transformed parameters {';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NSUB] id_terms;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NOCC] occ_terms;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NTRL] trl_terms;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NTID] trlxid_terms;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NOID] occxid_terms;');
            stan_in{end+1,1} = ...
                sprintf('  vector[NTO] trlxocc_terms;');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  id_terms = sig_id * id_raw;');
            stan_in{end+1,1} = ...
                sprintf('  occ_terms = sig_occ * occ_raw;');
            stan_in{end+1,1} = ...
                sprintf('  trl_terms = sig_trl * trl_raw;');
            stan_in{end+1,1} = ...
                sprintf('  trlxid_terms = sig_trlxid * trlxid_raw;');
            stan_in{end+1,1} = ...
                sprintf('  occxid_terms = sig_occxid * occxid_raw;');
            stan_in{end+1,1} = ...
                sprintf('  trlxocc_terms = sig_trlxocc * trlxocc_raw;');
        end
        
        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = '  ';
        stan_in{end+1,1} = 'model {';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  vector[NOBS] y_hat;');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  for(i in 1:NOBS)');
            stan_in{end+1,1} = ...
                sprintf('    y_hat[i] = pop_int + id_terms[id[i]] +');
            stan_in{end+1,1} = ...
                sprintf('    occ_terms[occ[i]] + trl_terms[trl[i]] +');
            stan_in{end+1,1} = ...
                sprintf('    trlxid_terms[trlxid[i]] + occxid_terms[occxid[i]] +');
            stan_in{end+1,1} = ...
                sprintf('    trlxocc_terms[trlxocc[i]];');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  meas ~ normal(y_hat,sig_err);');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  id_raw ~ normal(0,1);');
            stan_in{end+1,1} = ...
                sprintf('  occ_raw ~ normal(0,1);');
            stan_in{end+1,1} = ...
                sprintf('  trl_raw ~ normal(0,1);');
            stan_in{end+1,1} = ...
                sprintf('  trlxid_raw ~ normal(0,1);');
            stan_in{end+1,1} = ...
                sprintf('  occxid_raw ~ normal(0,1);');
            stan_in{end+1,1} = ...
                sprintf('  trlxocc_raw ~ normal(0,1);');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  sig_id ~ cauchy(0,20);');
            stan_in{end+1,1} = ...
                sprintf('  sig_occ ~ cauchy(0,.1);');
            stan_in{end+1,1} = ...
                sprintf('  sig_trl ~ cauchy(0,20);');
            stan_in{end+1,1} = ...
                sprintf('  sig_err ~ cauchy(0,20);');
        end
        
        stan_in{end+1,1} = '  ';
        
        for i=1:ndchunks
            stan_in{end+1,1} = ...
                sprintf('  sig_trlxid ~ cauchy(0,10);');
            stan_in{end+1,1} = ...
                sprintf('  sig_occxid ~ cauchy(0,10);');
            stan_in{end+1,1} = ...
                sprintf('  sig_trlxocc ~ cauchy(0,.05);');
        end
        
        stan_in{end+1,1} = '}';
        %store the cmdstan syntax
        REL.stan_in = stan_in;
        
        %create structure for the data to be sent to cmdstan
        data = struct;
        
        %prep data structure for cmdstan
        for i = 1:ndchunks
            fieldname = sprintf('NOBS');
            fieldvalue = height(darray.data{i});
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NOCC');
            fieldvalue = length(unique(darray.data{i}.time));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NSUB');
            fieldvalue = length(unique(darray.data{i}.id2));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NTRL');
            fieldvalue = length(unique(darray.data{i}.trl2));
            data.(fieldname) = fieldvalue;
        end
        
        
        for i = 1:ndchunks
            fieldname = sprintf('occ');
            fieldvalue = darray.data{i}.time2;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('id');
            fieldvalue = darray.data{i}.id2;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('trl');
            fieldvalue = darray.data{i}.trl2;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('meas');
            fieldvalue = darray.data{i}.meas;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('trlxid');
            fieldvalue = darray.data{i}.idxtrl;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('occxid');
            fieldvalue = darray.data{i}.idxtim;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('trlxocc');
            fieldvalue = darray.data{i}.trlxtim;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NTID');
            fieldvalue = length(unique(darray.data{i}.idxtrl));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NOID');
            fieldvalue = length(unique(darray.data{i}.idxtim));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NTO');
            fieldvalue = length(unique(darray.data{i}.trlxtim));
            data.(fieldname) = fieldvalue;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        
        
        %store the cmdstan syntax
        REL.stan_in = stan_in;
        
        %create structure for the data to be sent to cmdstan
        data = struct;
        
        %prep data structure for cmdstan
        for i = 1:ndchunks
            fieldname = sprintf('NOBS%d',i);
            fieldvalue = height(darray.data{i});
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NOCC%d',i);
            fieldvalue = ntime;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NSUB%d',i);
            fieldvalue = length(unique(id2));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NTRL%d',i);
            fieldvalue = length(unique(darray.data{i}.trl2));
            data.(fieldname) = fieldvalue;
        end
        
        
        for i = 1:ndchunks
            fieldname = sprintf('occ%d',i);
            fieldvalue = darray.data{i}.time2;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('id%d',i);
            fieldvalue = darray.data{i}.id2;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('trl%d',i);
            fieldvalue = darray.data{i}.trl2;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('meas%d',i);
            fieldvalue = darray.data{i}.meas;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('trlxid%d',i);
            fieldvalue = darray.data{i}.idxtrl;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('occxid%d',i);
            fieldvalue = darray.data{i}.idxtim;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('trlxocc%d',i);
            fieldvalue = darray.data{i}.trlxtim;
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NTID%d',i);
            fieldvalue = length(unique(darray.data{i}.idxtrl));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NOID%d',i);
            fieldvalue = length(unique(darray.data{i}.idxtim));
            data.(fieldname) = fieldvalue;
        end
        
        for i = 1:ndchunks
            fieldname = sprintf('NTO%d',i);
            fieldvalue = length(unique(darray.data{i}.trlxtim));
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
            'file_overwrite', true);
        
        
        %%%?????control = list(adapt_delta=.98, max_treedepth=15)
        
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

function convstats = era_storeconv(fit)
%pull out r_hats from Stanfit model

%pull summary using the print function
%there's no other way to get at the r_hat values, so output will also be
%printed in the command window
output = print(fit);

%define cell array that holds the r_hats
convstats = cell(0,3);
convstats{1,1} = 'name';
convstats{1,2} = 'n_eff';
convstats{1,3} = 'r_hat';

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

end
